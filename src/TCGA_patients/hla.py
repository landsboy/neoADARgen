"""
hla.py
------
Utility functions for reading and normalizing HLA alleles
from TCGA patient tables.
"""

from typing import List
import logging

log = logging.getLogger(__name__)


def find_HLA(patient_id: str, df) -> List[str]:
    """
    Extract and normalize all HLA alleles for a given TCGA patient
    from the HLA typing table.

    Parameters
    ----------
    patient_id : str
        TCGA patient ID (e.g., "TCGA-AB-1234").
    df : pandas.DataFrame
        DataFrame containing HLA typings for the current cancer type.
        The first column must contain patient IDs, and all other columns
        should contain HLA allele strings (e.g., "HLA-A*02:01", "HLA-B*44:03").

    Returns
    -------
    list of str
        List of cleaned HLA alleles in NetMHCpan-compatible format.
        Example: ["HLA-A0201", "HLA-B4403"]
    """
    # Find the row corresponding to the patient
    patient_row = df[df.iloc[:, 0] == patient_id]
    hla_list: List[str] = []

    if not patient_row.empty:
        # Take all values from column index 2 onward (HLA alleles)
        hlas = patient_row.iloc[:, 2:].values.flatten().tolist()

        for hla in hlas:
            if not isinstance(hla, str) or not hla.strip():
                continue
            parts = hla.split('*')
            if len(parts) == 2:
                # Convert e.g. "HLA-A*02:01" → "HLA-A0201"
                head = parts[0]
                tail = parts[1].replace(':', '')
                converted = f'{head}{tail}'
                hla_list.append(converted)

    if not hla_list:
        log.warning(f"No HLA alleles found for patient {patient_id}")

    return hla_list
