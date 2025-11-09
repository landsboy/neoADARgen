"""
preprocessing.py
----------------
Extracts the DNA sequence surrounding a genomic mutation,
adjusts the reading frame, and introduces the mutation into the sequence.

This module wraps external shell scripts (`find_seq_of_mutation.sh`, `find_position.sh`)
to locate the genomic region and frame information, then builds
the mutated DNA isoform for downstream NetMHCpan prediction.
"""

import os
import re
import logging
import subprocess
from typing import List, Tuple, Optional
from Bio.Seq import Seq
from .config import AppConfig

log = logging.getLogger(__name__)


# ============================================================
# Helper — pattern matching for mutation strings
# ============================================================
def checking_input(mutation_input: str) -> Optional[re.Match]:
    """
    Identify the type of mutation (substitution / deletion / insertion)
    and extract chromosomal components.

    Returns
    -------
    re.Match | None
        Regex match object containing:
        (chromosome, position, type/variant, sequence)
    """
    if "del" in mutation_input:
        pattern = r"(\w+):g\.(\d+)([a-z]+)([A-Z]+)"
    elif "ins" in mutation_input:
        pattern = r"(\w+):g\.(\d+)_\d+([a-z]+)([A-Z]+)"
    elif ">" in mutation_input:
        pattern = r"(\w+):g\.(\d+)(\w)>(\w)"
    else:
        return None
    return re.search(pattern, mutation_input)

def build_dna_seq(seq, mut_pos, row, mut_type, mut_seq, seq_len, strand, reading_frame, mutation_input, cfg) -> Optional[Tuple[str, str, List[int], int, int]]:
    """
    Build DNA sequence around a mutation while aligning it to the correct reading frame
    and respecting strand orientation and CDS boundaries.
    Returns the mutated DNA sequence, strand, position list, mutation length, and mutation position.
    """
    temp = seq[::-1]
    frame_init_pos = int(mut_pos) - int((seq_len/2) - 1)
    reading_frame = int(row[5])
    if strand == "+":
        cds_start_pos = int(row[1])
        cds_end_pos = int(row[2])
        count = 0
        count2 = 1
    else:
        count = 1
        count2 = 0
        seq = temp
        frame_init_pos = int(mut_pos) + int((seq_len/2) - 1)
        cds_start_pos = int(row[2])
        cds_end_pos = int(row[1])
    end = len(seq)
    flag = True
    if abs(int(mut_pos) - cds_start_pos) < int((seq_len/2) - 1): 
        start = int((seq_len/2) - 1) - abs(int(mut_pos) - cds_start_pos) + reading_frame + count
        frame_init_pos = cds_start_pos + reading_frame 
        flag = False
    if abs(cds_end_pos - int(mut_pos)) < int(seq_len/2):
        end = abs(cds_end_pos - frame_init_pos) 
        if not flag:
            end += start
    if flag:
        if (reading_frame == 0):
            seq_frame = abs(frame_init_pos - cds_start_pos + count) % 3 
        elif(reading_frame == 1):
            seq_frame = abs(frame_init_pos - 1 - cds_start_pos) % 3 
        else:
            seq_frame = abs(frame_init_pos - 2 - cds_start_pos - count) % 3

        if seq_frame == 0:
            start = seq_frame 
        elif seq_frame == 1:
            start = seq_frame + 1
        elif seq_frame == 2:
            start = seq_frame - 1
            
    codon = [seq[i:i+3] for i in range(start, end, 3)]

    # Creats mutatation in the seq for each isoform:
    mut_isoforms = []
    pos_list = []
    pos = int(seq_len/2) - start - count2 + count
    if strand == "-":
        if mut_type == "ins":
            mut_seq = mut_seq[::-1]
        pos = pos - len(mut_seq)

    # Convert the list of codons to a DNA sequence
    dna_seq = ''.join(codon)
    if mut_type == "del":
        # Remove the sequence from the DNA sequence:
        mut_dna_seq = dna_seq[:pos] + dna_seq[pos + len(mut_seq):]
        pos_list.append(pos)
    elif mut_type == "ins":
        pos = int(seq_len/2) - start 
        # Insert the sequence to the DNA sequence:
        mut_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos:]
        for i in range(len(mut_seq)):
            pos_list.append(pos + i)
    else:
        mut_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos + 1:]
        pos_list.append(pos)
    
    if str(Seq(mut_dna_seq).translate()) == str(Seq(dna_seq).translate()):  #TODO
        log.debug(f"The {mutation_input} is a synonymous mutation.")
        return None

    m = len(mut_seq)
    if mut_type == "del":
        m = 0

    if len(mut_dna_seq) < cfg.runtime.mer_length * 3:
        log.debug(f"Sequence too short for NetMHCpan: {mutation_input}")
        return None
    
    if strand == "-":
        mut_dna_seq = str(Seq(mut_dna_seq).complement())

    mut_isoforms.append(mut_dna_seq)



    return mut_isoforms, strand, pos_list, m, pos

# ============================================================
# Main Preprocessing Function
# ============================================================

def preprocessing(mutation: tuple, k: int, cfg: AppConfig) -> Optional[Tuple[List[str], str, List[int], int, int]]:
    """
    Preprocess a genomic mutation to extract its sequence context and isoform using pure Python tools.

    Parameters
    ----------
    mutation : tuple
        (mutation_string, gene_name, transcript_id, HGVSp_Short)
    k : int
        Mutation index (for temporary file naming)
    cfg : config object
        Configuration settings (must include runtime.mer_length and runtime.num_nuc_around_mut)

    Returns
    -------
    (mut_isoforms, strand, pos_list, mut_len, pos) or None
    """

    # ==========================
    # Parse mutation
    # ==========================
    mutation_input = mutation[0]
    match = checking_input(mutation_input)
    if not match:
        log.warning(f"Skipping unrecognized mutation format: {mutation_input}")
        return None

    mut_chr = match.group(1)
    mut_pos = int(match.group(2))
    mut_type = "Substitution" if ">" in mutation_input else match.group(3)
    mut_seq = match.group(4)

    # ==========================
    # Extract sequence context
    # ==========================
    window = cfg.runtime.num_nuc_around_mut
    seq_len = 2 * window + 2 * 26 + 1   
    fasta_path = cfg.paths.hg38_fa
    bed_path = os.path.join(cfg.paths.sup_dir, "bed6_of_ensemble.bed")

    cmd = (
        f'sh {os.path.join("src", "find_position.sh")} '
        f'{mut_chr} {mut_pos} {k} {bed_path}'
    )

    transcript_proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    res_list = [line for line in transcript_proc.stdout.split("\n") if line.strip()]
    if not res_list:
        log.error(f"No transcript info found for {mutation_input}")
        return None
    
    row = next(
        (r.split("\t") for r in res_list if mutation[2] in r.split("\t")[3]),
        res_list[0].split("\t")
    )

    strand = None
    reading_frame = None
    # Extract STRAND information (+ / -) Usually stored in column 5 (index 4)
    if len(row) >= 5 and row[4] in {"+", "-"}:
        strand = row[4]

    # Extract READING FRAME Often stored in column 6 (index 5) or column 7 (index 6)
    if len(row) >= 6 and row[5].isdigit():
        reading_frame = int(row[5])
    elif len(row) >= 7 and row[6].isdigit():
        reading_frame = int(row[6])

    # Error handling for missing strand
    if strand is None:
        log.error(
            f"Missing or invalid strand for {mutation[0]} | row={row}"
        )
        return None

    # Error handling for missing reading frame
    if reading_frame is None:
        log.error(
            f"Missing or invalid reading frame (RF) for {mutation[0]} "
            f"(transcript {mutation[2]}). Row fields: {row}"
        )
        return None

    half = seq_len // 2
    if strand == "+":
        start_seq = mut_pos - half
        end_seq = mut_pos + half + 1
    else:
        start_seq = mut_pos - half - 1
        end_seq = mut_pos + half
    
    cmd = (
        f'sh {os.path.join(os.getcwd(), "src", "find_seq_of_mutation.sh")} '
        f'{mut_chr} {start_seq} {end_seq} {k} {fasta_path}'
    )
    seq_proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    seq = seq_proc.stdout.replace("\n", "")
    if not seq:
        log.error(f"Failed to extract sequence for {mutation_input}")
        return None

    return build_dna_seq(seq, mut_pos, row, mut_type, mut_seq, seq_len, strand, reading_frame, mutation_input, cfg)