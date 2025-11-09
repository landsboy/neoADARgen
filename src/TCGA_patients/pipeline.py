"""
pipeline.py
------------
Main workflow controller for the neoADARgen pipeline.
Handles mutation parsing, ADAR editing simulation,
NetMHCpan predictions, and output file generation.
"""

import os
import logging
from typing import List, Dict
import pandas as pd
from Bio.Seq import Seq

# Internal imports
from .config import AppConfig
from .preprocessing import preprocessing
from .mutations import create_list_of_mutations, generate_single_edit_sequences, generate_double_edit_sequences
from .mhcpan import run_netmhcpan, parse_netmhcpan_output
from .hla import find_HLA
from .io_utils import find_duplicate_patients_in_data, add_TPM

# Initialize logger
log = logging.getLogger(__name__)


# =============================================================================
# MAIN PIPELINE
# =============================================================================
def main(cfg: AppConfig) -> None:
    """
    Main entry point for the neoADARgen pipeline.

    Iterates over all cancer types in the project directory,
    loads HLA typing data, processes patients, and writes final results.
    """
    project_dir = cfg.paths.project_dir
    results_dir = cfg.paths.results_dir
    sup_dir = cfg.paths.sup_dir

    os.makedirs(results_dir, exist_ok=True)
    log.info(f"Starting neoADARgen pipeline in: {project_dir}")

    for primary_site in os.listdir(project_dir):
        site_dir = os.path.join(project_dir, primary_site)
        if not os.path.isdir(site_dir):
            continue

        os.makedirs(os.path.join(results_dir, primary_site), exist_ok=True)
        hla_path = os.path.join(sup_dir, "HLA", f"{primary_site}.tsv")

        if not os.path.exists(hla_path):
            log.warning(f"HLA file missing for cancer type: {primary_site}")
            continue

        hla_df = pd.read_csv(hla_path, sep="\t")
        log.info(f"Processing cancer type: {primary_site}")

        duplicates, patient_dict = find_duplicate_patients_in_data(primary_site, project_dir)
        patient_folders = [
            f for f in os.listdir(site_dir)
            if os.path.isdir(os.path.join(site_dir, f))
        ]

        log.debug(f"Found {len(patient_folders)} patients for {primary_site}")
        run_patient_batch(patient_folders, cfg, primary_site, duplicates, hla_df, patient_dict)

        # Add TPM data after processing all patients
        add_TPM(primary_site, duplicates, cfg, results_dir)

    log.info("neoADARgen pipeline completed successfully.")


# =============================================================================
# PATIENT BATCH PROCESSING
# =============================================================================
def run_patient_batch(
    patient_folders: List[str],
    cfg: AppConfig,
    primary_site: str,
    duplicates: List[str],
    hla_df: pd.DataFrame,
    patient_dict: Dict[str, str],
) -> None:
    """
    Process all patients for a specific cancer type.
    Each patient folder is handled independently.
    """
    for idx, folder in enumerate(patient_folders):
        try:
            mut_list, full_patient_id, n_mutations = create_list_of_mutations(
                os.path.join(cfg.paths.project_dir, primary_site, folder)
            )
        except Exception as e:
            log.error(f"Failed to read mutations for {folder}: {e}", exc_info=True)
            continue

        patient_id = "-".join(full_patient_id.split("-")[0:3])
        res_path = define_output_path(cfg, primary_site, patient_id, full_patient_id, duplicates, patient_dict)
        hla_list = list(set(find_HLA(patient_id, hla_df)))

        if not hla_list:
            log.warning(f"No HLA found for patient {patient_id} in {primary_site} project, skipping.")
            continue

        log.info(f"Processing patient {patient_id} with {len(mut_list)} CDS mutations.")

        with open(res_path, "w") as fout:
            write_file_header(fout, patient_id, folder, n_mutations, len(mut_list))

            for mut in mut_list:
                process_mutation(fout, mut, idx, hla_list, cfg)


# ============================================================
#  Process single mutation
# ============================================================
def process_mutation(fout, mutation, idx, hla_list, cfg: AppConfig):
    """Master function to process a single mutation and write best neoantigen result."""
    result = preprocessing(mutation, idx, cfg)
    if result is None:
        return

    mut_isoforms, strand, pos_list, mut_len, pos = result

    for edit_option in cfg.runtime.edit_modes:
        without_editing, single_editing, double_editing = edit_option_initialization(edit_option)
        all_results = []

        if without_editing:
            all_results = run_mhc_binding_analysis(mut_isoforms, cfg, hla_list, idx)
            seq = list(mut_isoforms[0])
            gRNA = "".join(seq[max(pos - cfg.runtime.num_nuc_around_mut, 0):
                            min(pos + cfg.runtime.num_nuc_around_mut + mut_len, len(seq))])
            write_best_result(fout, mutation, all_results, mut_isoforms, strand, pos, mut_len, edit_option, cfg, gRNA, None)

        if single_editing or double_editing:
            adar_seqs = []
            if single_editing:
                adar_seqs += generate_single_edit_sequences(mut_isoforms, pos, mut_len, pos_list, cfg)
            if double_editing:
                adar_seqs += generate_double_edit_sequences(mut_isoforms, pos, mut_len, pos_list, cfg)
            adar_unique = list(set(adar_seqs))
            all_results += run_mhc_binding_analysis([a[0] for a in adar_unique], cfg, hla_list, idx)
            write_best_result(fout, mutation, all_results, mut_isoforms, strand, pos, mut_len, edit_option, cfg, None, adar_unique)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def write_best_result(fout, mutation, sb_wb_results, mut_isoforms, strand, pos, m, edit_option, cfg, gRNA , unique_ADAR_mut_list):
    """
    Write the strongest predicted neoantigen for the current mutation
    to the output results file.
    """
    if not sb_wb_results:
        seq_list = list(mut_isoforms[0])
        gRNA_NA = "".join(
            seq_list[max(pos - cfg.runtime.num_nuc_around_mut, 0):min(pos + cfg.runtime.num_nuc_around_mut + m, len(seq_list))]
        )
        fout.write(
            f"{mutation[1]}\t{mutation[2]}\t{mutation[0]}\t{strand}\t{mutation[3]}"
            f"\tNA\t{gRNA_NA}\tNA\tNA\tNA\tNA\t{edit_option}\n"
        )
        return

    best = sorted(sb_wb_results, key=lambda x: (x[2], x[4], x[5]))[0]
    isoform_index = best[0]
    best_peptide = best[1]

    # Select DNA sequence depending on edit mode
    if edit_option in [1, 2]:
        sb_seq = Seq(unique_ADAR_mut_list[isoform_index][0])
        guide_name = unique_ADAR_mut_list[isoform_index][1]
    else:
        # sb_seq = Seq(NEW_mut_isoforms[isoform_index])
        sb_seq = Seq(mut_isoforms[isoform_index])
        guide_name = gRNA

    # Translate to protein
    protein = sb_seq.translate()
    protein_str = str(protein)

    # Replace stop codons (*) with X
    if "*" in protein_str:
        protein_str = protein_str.replace("*", "X")
        protein = Seq(protein_str)

    # Find peptide in protein
    pos_p = protein.find(best_peptide)
    if pos_p == -1:
        ASO_target = "NA"
    else:
        ASO_target = sb_seq[(pos_p * 3):((pos_p + len(best_peptide)) * 3)]

    # Pack final result tuple 
    best_ASO = (
        guide_name,          # gRNA / ADAR guide
        ASO_target,          # DNA target sequence
        best[2],     # rank
        mutation,            # full mutation info
        best[3],     # HLA
        strand,              # +/-
        best[4],     # rank_ba
        best[5],     # affinity
    )

    # Write output line 
    fout.write(
        f"{best_ASO[3][1]}\t{best_ASO[3][2]}\t{best_ASO[3][0]}\t{best_ASO[5]}\t"
        f"{best_ASO[3][3]}\t{str(best_ASO[1])}\t{str(best_ASO[0])}\t{best_ASO[2]}\t"
        f"{best_ASO[6]}\t{best_ASO[7]}\t{best_ASO[4]}\t{edit_option}\n"
    )


def define_output_path(cfg, primary_site, patient_id, file_id, duplicates, patient_dict) -> str:
    """
    Define the output path for a patient results file,
    taking into account duplicated patient samples.
    """
    if patient_id in duplicates:
        if patient_id not in patient_dict or file_id != patient_dict[patient_id]:
            return ""
        return os.path.join(cfg.paths.results_dir, primary_site, f"{file_id}.tsv")
    return os.path.join(cfg.paths.results_dir, primary_site, f"{patient_id}.tsv")


def write_file_header(fout, patient_id, file_id, n_mutations, n_cds):
    """
    Write a formatted header to the per-patient results file.
    """
    fout.write("# neoADARgen 1.0\n")
    fout.write(f"# Patient ID: {patient_id}\n")
    fout.write(f"# File ID: {file_id}\n")
    fout.write(f"# Num of mutations: {n_mutations}\n")
    fout.write(f"# Num of CDS mutations: {n_cds}\n")
    fout.write(
        "Gene_Name\tTranscript_ID\tMutation\tStrand\tHGVSp_Short\tBest-target\tGuide-RNA\t"
        "Rank(%)\tAffinity(nM)\tRank_BA(%)\tHLA\tEDITS\n"
    )


def edit_option_initialization(edit_option: int) -> tuple[bool, bool, bool]:
    """
    Initialize the editing mode flags for a given edit option.
    """
    if edit_option == 0:
        return True, False, False      # No ADAR editing
    elif edit_option == 1:
        return False, True, False      # Single A→G edit
    else:
        return False, False, True      # Double A→G edits


def create_fasta_file(temp_path, list_of_seq, k):
    """
    Create a FASTA file from a list of sequences for NetMHCpan input.
    """
    with open(os.path.join(temp_path, "TEMP", f"input_seq{k}.fsa"), "w") as fasta_file:
        for i, sequence in enumerate(list_of_seq):
            fasta_file.write(f">{i}\n{str(Seq(sequence).translate())}\n")


def run_mhc_binding_analysis(seqs: List[str], cfg, hla_list, idx):
    """Run NetMHCpan on edited sequences and collect predictions."""
    sb_wb_tuple = []
    if seqs:
        create_fasta_file(cfg.paths.sup_dir, seqs, idx)
        for hla in hla_list:
            output = run_netmhcpan(cfg, hla, idx)
            sb_wb_tuple.extend(parse_netmhcpan_output(output))
    return sb_wb_tuple