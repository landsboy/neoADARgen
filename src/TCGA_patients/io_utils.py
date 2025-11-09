"""
io_utils.py
-----------
Utility functions for file I/O operations used across TCGA_patients.
Includes helpers for reading gzip-compressed files and scanning directories.
"""
import shutil
import os
import csv
from pathlib import Path
import gzip
from typing import Optional, Dict, List, Tuple
import logging
from .config import AppConfig

log = logging.getLogger(__name__)


def setup_environment(netmhc_path: Path):
    """
    Configure environment variables required for NetMHCpan execution.
    Uses the 'netmhc_path' from config instead of hard-coded paths.
    """
    os.environ["PLATFORM"] = "Linux_x86_64"
    os.environ["NMHOME"] = str(netmhc_path)
    os.environ["NETMHCpan"] = str(netmhc_path / "Linux_x86_64")  # or adjust if different
    os.environ["NetMHCpanWWWPATH"] = "/services/NetMHCpan/tmp/"
    os.environ["NetMHCpanWWWDIR"] = "/usr/opt/www/pub/CBS/services/NetMHCpan/tmp"
    os.environ["TMPDIR"] = "/tmp"
    os.environ["DTUIBSWWW"] = "www"

    print("[INFO] Environment variables for NetMHCpan were successfully configured.")

def first_gz_in_dir(dir_path: Path) -> Optional[Path]:
    """
    Return the first .gz file found in a given directory.

    Parameters
    ----------
    dir_path : Path
        Directory path to search.

    Returns
    -------
    Path or None
        Path to the first .gz file found, or None if none exist.
    """
    for p in dir_path.iterdir():
        if p.suffix == ".gz" and p.is_file():
            log.debug(f"Found GZ file: {p}")
            return p
    log.warning(f"No .gz files found in {dir_path}")
    return None

def find_duplicate_patients_in_data(primary_site: str, project_path: Path) -> Tuple[List[str], Dict[str, str]]:
    """
    Identify TCGA patients that have multiple mutation files (duplicate entries).

    Parameters
    ----------
    primary_site : str
        TCGA cancer project code (e.g. BRCA, SKCM, OV).
    project_path : Path
        Path to the root TCGA project directory.

    Returns
    -------
    dup_pat : list of str
        List of patient IDs that have two or more mutation files.
    new_dic : dict
        Mapping of patient_id → selected file_id (one file per patient).
        The selection is based on tumor type and UUID mapping rules.
    """

    all_patient: Dict[str, List[str]] = {}

    # Iterate through patient folders
    for file_id in [
        f for f in os.listdir(os.path.join(project_path, primary_site))
        if os.path.isdir(os.path.join(project_path, primary_site, f))
    ]:
        patient_dir = os.path.join(project_path, primary_site, file_id)

        # Find the first .gz mutation file
        gz_file_path = None
        for file_name in os.listdir(patient_dir):
            if file_name.endswith(".gz"):
                gz_file_path = os.path.join(patient_dir, file_name)
                break

        if not gz_file_path:
            log.warning(f"No .gz file found for {file_id} in {patient_dir}")
            continue

        # Extract patient ID from the first mutation entry
        try:
            with gzip.open(gz_file_path, "rt") as file:
                for idx, line in enumerate(file):
                    if idx <= 8:
                        continue
                    columns = line.strip().split("\t")
                    patient_key = "-".join(columns[15].split("-")[0:3])
                    all_patient.setdefault(patient_key, []).append(columns[15])
                    break
        except Exception as e:
            log.error(f"Failed to read {gz_file_path}: {e}")
            continue

    # Identify duplicates
    dup_pat = [key for key, value in all_patient.items() if len(value) >= 2]

    # Read the UUID-to-CASE reference file
    uuid_file = os.path.join(project_path, primary_site, "UUIDtoCASE.tsv")
    if not os.path.exists(uuid_file):
        log.warning(f"UUIDtoCASE.tsv not found in {primary_site}")
        return dup_pat, {}

    with open(uuid_file, "r") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        data = [row[0] for row in reader]

    # Resolve duplicates by tumor sample preference
    new_dic: Dict[str, str] = {}
    for key, val in all_patient.items():
        if key in dup_pat:
            for v in val:
                if (
                    v.split("-")[3][:2] == "01"
                    and v.split("-")[4][-1] == "D"
                    and "-".join(v.split("-")[0:4]) in data
                ):
                    new_dic[key] = v
            if key not in new_dic:
                for v in val:
                    if (
                        v.split("-")[3][:2] == "01"
                        and "-".join(v.split("-")[0:4]) in data
                    ):
                        new_dic[key] = v

    log.info(f"{len(dup_pat)} duplicate patients found in {primary_site}")
    return dup_pat, new_dic

def add_TPM(primary_site: str, duplicate_patients: list[str], cfg: AppConfig, res_path: Path) -> None:
    """
    Append TPM expression values to each patient result file.

    Parameters
    ----------
    primary_site : str
        TCGA project code (e.g., BRCA, SKCM, OV).
    duplicate_patients : list of str
        Patients that have multiple mutation files.
    project_path : Path
        Root path of the project (contains mutation folders).
    res_path : Path
        Results output directory.
    """

    uuid_path = Path(cfg.paths.project_dir) / primary_site / "UUIDtoCASE.tsv"
    if not uuid_path.exists():
        log.warning(f"UUIDtoCASE.tsv not found for {primary_site}")
        return

    # Read the UUID mapping
    with open(uuid_path, "r") as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        data = [(row[1], row[2], row[0]) for row in reader] # (Case ID, TPM folder name, BCR sample UUID)

    for filename in os.listdir(Path(res_path) / primary_site):
        mapping = {}

        # Extract patient ID correctly (handles cases with multiple hyphens)
        base_name = filename.split(".")[0]
        if base_name.count("-") > 3:
            parts = base_name.split("-")
            patient_id = "-".join(parts[:4])
            temp = "-".join(parts[:3])
        else:
            patient_id = base_name
            temp = base_name

        # Match patient to UUID mapping row
        for item in data:
            if (patient_id not in duplicate_patients and item[0] == patient_id) or (
                temp in duplicate_patients and item[2] == patient_id):
                
                # Path to TPM folder for this sample
                folder_path = Path(cfg.paths.sup_dir) / "TPM" / primary_site / item[1]

                if not folder_path.exists():
                    log.warning(f"TPM folder not found: {folder_path}")
                    break

                # Load TPM file
                tsv_files = [f for f in os.listdir(folder_path) if f.endswith(".tsv")]
                if not tsv_files:
                    log.warning(f"No TPM file in {folder_path}")
                    break

                tpm_file_path = folder_path / tsv_files[0]
                with open(tpm_file_path, "r") as tpm_file:
                    for i, line in enumerate(tpm_file):
                        if i >= 6:  # skip header lines
                            columns = line.strip().split("\t")
                            key = columns[6]   # Gene name
                            value = columns[1] # TPM value
                            mapping[value] = key

                input_file = Path(res_path) / primary_site / filename
                temp_file = input_file.with_suffix(".tmp")

                with open(input_file, "r") as patient_file, open(temp_file, "w") as new_patient_file:
                    for i, line in enumerate(patient_file):
                        if i >= 6:
                            columns = line.strip().split("\t")
                            gene = columns[0]
                            tpm_value = mapping.get(gene, "")
                            if tpm_value:
                                new_patient_file.write(line.strip() + f"\t{tpm_value}\n")
                            else:
                                new_patient_file.write(line)
                        elif i == 5:
                            new_patient_file.write(line.strip() + "\tTPM\n")
                        else:
                            new_patient_file.write(line)
                
                # Replace original file atomically
                shutil.move(temp_file, input_file)            
                break

    log.info(f"TPM values added for {primary_site}")