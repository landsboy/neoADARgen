"""
mhcpan.py
---------
Functions for running and parsing NetMHCpan output.
This module provides wrappers for executing the NetMHCpan binary
and extracting peptide-binding results in a structured format.
"""

import subprocess
from typing import List, Tuple
import logging
from .config import AppConfig
import os

log = logging.getLogger(__name__)


def run_netmhcpan(cfg: AppConfig, hla: str, idx: int) -> str:
    """
    Run the NetMHCpan tool on a given FASTA file using configuration settings.

    Parameters
    ----------
    cfg : AppConfig
        Configuration object containing paths and runtime parameters.
    hla : str
        HLA allele (e.g., "HLA-A02:01") to use in NetMHCpan prediction.
    idx : int
        Unique index for naming temporary files.

    Returns
    -------
    str
        Raw text output from NetMHCpan.
    """
    bin_path = cfg.paths.netmhc_path / "Linux_x86_64" / "bin" / "netMHCpan"
    cmd = f"{bin_path} -f {os.path.join(cfg.paths.sup_dir, 'TEMP', f'input_seq{idx}.fsa')} -a {hla} -l {cfg.runtime.mer_length} -BA -t 2.1"
    log.debug(f"Running NetMHCpan command: {cmd}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"NetMHCpan failed: {result.stderr}")
        return ""
    return result.stdout


def parse_netmhcpan_output(text: str) -> List[Tuple[int, str, str, str, str, str]]:
    """
    Parse NetMHCpan output text and extract relevant peptide-binding information.

    Returns
    -------
    list of tuples:
        (index, peptide, rank, hla, rank_ba, affinity)
    """
    results: List[Tuple[int, str, str, str, str, str]] = []

    for line in text.splitlines():
        if not line or line.startswith("#") or line.startswith("Pos"):
            continue
        tokens = [t for t in line.split(" ") if t.strip()]
        if len(tokens) < 16:
            continue
        flag = tokens[-1]
        if flag not in {"SB", "WB"}:  # SB = strong binder, WB = weak binder
            continue

        try:
            idx = int(tokens[10])
        except ValueError:
            continue

        peptide = tokens[2]
        hla = tokens[1].replace("*", "")
        rank = tokens[12]
        rank_ba = tokens[15]
        affinity = tokens[14]

        results.append((idx, peptide, rank, hla, rank_ba, affinity))

    return results


