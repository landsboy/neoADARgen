# neoADARgen
![Python](https://img.shields.io/badge/Python-3.10-blue)
![Shell](https://img.shields.io/badge/Shell-Bash-lightgrey)
![Conda](https://img.shields.io/badge/Conda-environment-44A833)
![NetMHCpan](https://img.shields.io/badge/NetMHCpan-4.1-FF6C37)
![Genome](https://img.shields.io/badge/Genome-hg38-2E8B57)
![RNAEditing](https://img.shields.io/badge/RNA_Editing-ADAR_A→I-8A2BE2)
![TCGA](https://img.shields.io/badge/Data-TCGA-4B0082)

A computational pipeline for generating neo-antigens through RNA editing.


neoADARgen is a bioinformatics tool designed to identify and engineer personalized tumor-specific neoantigens (editopes) by simulating A-to-I RNA editing events on somatic mutations.
The tool integrates mutation annotation, sequence extraction, RNA editing simulation, peptide generation, and MHC binding prediction via NetMHCpan 4.1.


This demo version of the tool in fact ran on three specific projects in TCGA (BRCA, GBN, SCKM), if you want to run on other projects you must download them from [TCGA](https://portal.gdc.cancer.gov/analysis_page?app=Downloads) and put them in the testdata folder.

An accompanying repo to the paper:

****Employing RNA editing to engineer personalized tumor-specific neoantigens (editopes)****

# Overview
The pipeline performs the following steps:

 - **Mutation parsing:** Reads patient-specific mutation annotations (MAF format).

 - **Sequence extraction:** Retrieves reference DNA sequences (±20 bp around each mutation).

 - **RNA editing simulation:** Applies A-to-I (read as G) edits in single or double combinations.

 - **Peptide generation:** Translates mutated and edited sequences into peptides of defined lengths (default: 9-mer).

 - **MHC binding prediction:** Uses NetMHCpan 4.1 to predict peptide–HLA binding affinities.

 - **Result summarization:** Outputs ranked neoantigen candidates per mutation.

# Setup
Requires locally installed version of [NetMHCpan4.1](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCpan&version=4.1&packageversion=4.1b&platform=Linux)

And requires a local download of the [human genome v.38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
 
1. Clone the repository
```
git clone https://github.com/landsboy/neo-ADARtigen.git
cd neoADARgen
```

2. Create conda environment
```
conda env create -f TCGA_environment.yml
conda activate TCGA_patients_env 
```
3. Running the Pipeline
The easiest way to run neoADARgen is by providing a configuration file (.yml) that defines all required paths and runtime parameters.
for example:
```yaml
paths:
  project_dir: "testdata"          # Directory containing raw patient mutation data
  results_dir: "results"           # Directory where all pipeline outputs will be saved
  sup_dir: "sup"                   # Directory with supplementary annotation files (e.g. HLA, TPM)
  netmhc_path: <path/to/your/netMHCpan4.1>
  hg38_fa: <path/to/your/hg38.fa>
  
runtime:
  edit_modes: [0, 1, 2]            # 0 = no RNA editing, 1 = single A→G editing, 2 = double editing
  mer_length: 9                    # Peptide length for NetMHCpan prediction (9-mer is the default)
  num_nuc_around_mut: 20           # Number of nucleotides to extract on each side of the mutation
  verbose: false                   # If true, enable detailed DEBUG logging
  log_file: "logs/TCGA_patients.log"   # Name of the log file saved in the results directory
```

Once the configuration file is ready (e.g. TCGA_config.yml), you can run the full pipeline with a single command:
```
python -m src.TCGA_patients.cli -c TCGA_config.yml
```
If you prefer, you can skip the config file and pass the parameters directly via the command line:
```
python -m src.TCGA_patients.cli \
  --project_dir testdata \
  --results_dir results \
  --sup_dir sup \
  --netmhc_path /path/to/netMHCpan-4.1 \
  --hg38_fa /path/to/hg38.fa \
  --verbose
  ```
  
  # Example Output

  For each patient, the pipeline generates an individual results file named according to their patient ID (e.g.results/BRCA/TCGA-AC-A2FK.tsv).
  
  ![Example output](sup/example.png)

  Within each file, all somatic mutations located in coding regions (CDS) are analyzed under three distinct conditions:

  0. Without RNA editing (original tumor mutation)

  1. A single A→G editing event (simulating ADAR activity at one site)

  2. With double A→G editing events (simulating ADAR activity at two positions)

  Each combination is processed through the NetMHCpan predictor to evaluate its HLA-binding affinity and neoantigen potential.

  This allows quantifying, for every patient, how RNA editing may increase the likelihood of generating strong-binding neoantigens — revealing novel tumor-specific “editopes”.

  # Getting help
  If you need help of any kind, feel free to open a new issue.