"""
config.py
----------
Defines the configuration schema for the TCGA_patients pipeline using Pydantic.
"""
from pathlib import Path
from typing import List
from pydantic import BaseModel, Field
import yaml

class Paths(BaseModel):
    project_dir: Path
    results_dir: Path
    sup_dir: Path
    netmhc_path: Path
    hg38_fa: Path

class Runtime(BaseModel):
    edit_modes: List[int] = Field(default_factory=lambda: [0, 1, 2])
    mer_length: int = 9
    num_nuc_around_mut: int = 20
    verbose: bool = False
    log_file: str = "TCGA_patients.log"

    @property
    def size_of_seq(self):
        return (2 * self.num_nuc_around_mut) + (2 * 26) + 1


class AppConfig(BaseModel):
    paths: Paths
    runtime: Runtime = Runtime()

def load_config(path: Path) -> AppConfig:
    with open(path, "r") as f:
        data = yaml.safe_load(f)
    cfg = AppConfig(**data)
    cfg.paths.results_dir.mkdir(parents=True, exist_ok=True)
    return cfg
