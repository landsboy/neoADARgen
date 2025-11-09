"""
mutations.py
------------
This module handles all mutation-level logic for the TCGA_patients pipeline:
- Parsing MAF files
- Constructing peptide windows
- Generating ADAR editing variants (single/double)
- Running NetMHCpan for MHC binding predictions
- Selecting and writing the best neoantigen candidate
"""

import gzip
import logging
from pathlib import Path
from typing import List, Tuple
from itertools import combinations
from .io_utils import first_gz_in_dir


log = logging.getLogger(__name__)


# ============================================================
#  Exception
# ============================================================
class MutationFileError(Exception):
    """Raised when a MAF file cannot be read or parsed properly."""
    pass

# ============================================================
#  Parse mutation list from MAF (.gz)
# ============================================================

def create_list_of_mutations(mut_dir: str | Path) -> Tuple[List[tuple], str, int]:
    """
    Extract mutation records from a patient's MAF .gz file.

    Parameters
    ----------
    mut_dir : str | Path
        Directory of a single patient that contains a MAF .gz file.

    Returns
    -------
    list_of_mutations : list of tuples
        (mutation_string, gene_symbol, transcript_id, HGVSp_Short)
    file_id : str
        Patient file/sample ID (from columns[15]).
    n_mutations : int
        Total number of mutation lines (after skipping headers).

    Raises
    ------
    MutationFileError
        If the MAF file is missing, unreadable, or malformed.
    """
    mut_dir = Path(mut_dir)
    gz_path = first_gz_in_dir(mut_dir)

    if gz_path is None:
        raise MutationFileError(f"No .gz MAF file found in {mut_dir}")

    list_mutations: List[tuple] = []
    line_count = 0
    file_id = ""

    try:
        with gzip.open(gz_path, "rt") as fh:
            for line in fh:
                line_count += 1
                # Start extracting mutations from the 7th line
                if line_count <= 8:
                    continue

                columns = line.rstrip("\n").split("\t")

                # Expecting at least 38 columns in the MAF file
                if len(columns) < 38:
                    raise MutationFileError(f"Malformed MAF line in {gz_path}: insufficient columns ({len(columns)})")

                file_id = columns[15]
                if not file_id:
                    raise MutationFileError(f"Missing file ID in {gz_path}")

                # Skip non-CDS mutations
                if not columns[35]:
                    continue

                # Skip mitochondrial chromosome
                if columns[4] == "chrM":
                    continue

                hgvsp_short = columns[36].split(".")[-1]
                if "=" in hgvsp_short:
                    continue

                # Skip nonsense and problematic frameshifts
                if "*" in hgvsp_short:
                    if "fs" in hgvsp_short:
                        num_aa = hgvsp_short.split("*")[-1]
                        if num_aa == "?" or (num_aa.isdigit() and int(num_aa) < 7):
                            continue
                        if not num_aa.isdigit():
                            continue
                    else:
                        continue

                # Build mutation string
                chrom = columns[4]
                start = columns[5]
                end = columns[6]
                var_type = columns[9]
                ref = columns[10]
                alt = columns[12]

                if var_type == "SNP":
                    new_mutation = f"{chrom}:g.{start}{ref}>{alt}"
                elif var_type == "DEL":
                    new_mutation = f"{chrom}:g.{start}del{ref}"
                elif var_type == "INS":
                    new_mutation = f"{chrom}:g.{start}_{end}ins{alt}"
                else:
                    continue

                list_mutations.append((new_mutation, columns[0], columns[37], hgvsp_short))

        n_mutations = max(line_count - 8, 0)

        if not list_mutations:
            raise MutationFileError(f"No valid mutations found in {gz_path}")

        return list_mutations, file_id, n_mutations

    except OSError as e:
        # Problem with reading gzip or file not found
        raise MutationFileError(f"Error reading {gz_path}: {e}") from e
    except Exception as e:
        # Any other parsing issue
        raise MutationFileError(f"Unexpected error parsing {gz_path}: {e}") from e

# ============================================================
#  Single-edit mode
# ============================================================
def generate_single_edit_sequences(mut_isoforms, pos, m, pos_list, cfg):
    """Generate all possible single A>G edits (ADAR)."""
    adar_list = []
    for isoform in mut_isoforms:
        for i in range(max(pos - cfg.runtime.num_nuc_around_mut + 1, 1),
                       min(pos + cfg.runtime.num_nuc_around_mut + m - 1, len(isoform) - 1)):
            seq = list(isoform)
            if seq[i] != "A" or i in pos_list:
                continue
            new_seq = seq.copy()
            new_seq[i] = "G"
            temp = isoform[:i] + "*" + isoform[i + 1:]
            gRNA = temp[max(pos - cfg.runtime.num_nuc_around_mut, 0):
                        min(pos + cfg.runtime.num_nuc_around_mut + m, len(seq))]
            if i % 3 == 0:
                start_index = max(i - 24, 0)
                end_index = i + 27
            elif i % 3 == 1:
                start_index = max(i - 25, 0)
                end_index = i + 26
            else:
                start_index = max(i - 26, 0)
                end_index = i + 25
            adar_list.append(("".join(new_seq[start_index:end_index]), gRNA, isoform[start_index:end_index]))
    return adar_list


# ============================================================
#  Double-edit mode
# ============================================================
def generate_double_edit_sequences(mut_isoforms, pos, m, pos_list, cfg):
    """Generate all possible double A>G edits (ADAR) within MER_LENGTH distance."""
    adar_pairs = []
    for isoform in mut_isoforms:
        seq = list(isoform)
        edit_sites = [i for i in range(max(pos - cfg.runtime.num_nuc_around_mut + 1, 1),
                                       min(pos + cfg.runtime.num_nuc_around_mut + m - 1, len(seq) - 1))
                      if seq[i] == "A" and i not in pos_list]

        for p1, p2 in combinations(edit_sites, 2):
            a = p2 - p1 + p1 % 3 + 1
            # The distance between two editing sites in amino acids
            AA_distance = -(-a // 3)
            # We want to produce a n-mer that has both editing sites
            if AA_distance > cfg.runtime.mer_length:
                continue
            new_seq = seq.copy()
            new_seq[p1] = "G"
            new_seq[p2] = "G"
            temp = "".join(new_seq)
            x = cfg.runtime.mer_length - AA_distance
            groups = [temp[i:i + 3] for i in range(0, len(temp), 3)]
            start_index = max(-(- (p1 + 1) // 3) - 1 - x, 0)
            end_index = -(- (p2 + 1) // 3) + x
            new_seq[p1] = "*"
            new_seq[p2] = "*"
            temp = "".join(new_seq)
            gRNA = temp[max(pos - cfg.runtime.num_nuc_around_mut, 0):
                        min(pos + cfg.runtime.num_nuc_around_mut + m, len(seq))]
            adar_pairs.append(("".join(groups[start_index:end_index]), gRNA, isoform[start_index * 3:end_index * 3]))
    return adar_pairs
