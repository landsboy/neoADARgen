# neoADARgen
A tool for creating a neo-antigen by RNA editing.

This demo version of the tool in fact ran on three specific projects in TCGA (BRCA, GBN, SCKM), if you want to run on other projects you must download them from [TCGA](https://portal.gdc.cancer.gov/analysis_page?app=Downloads) and put them in the testdata folder.

An accompanying repo to the paper:

****NAME OF THE ARTICLE****

# Getting help
If you need help of any kind, feel free to open a new issue.

# Setup
Requires locally installed version of [NetMHCpan4.1](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCpan&version=4.1&packageversion=4.1b&platform=Linux)

And requires a local download of the [human genome v.38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
 
1. Clone the repository
```
git clone https://github.com/landsboy/neo-ADARtigen.git
```

2. Create conda environment
```
conda env create -f neoADARgen.yml
conda activate neoADARgen  # test successful creation
```
3. Run the following file with the following arguments
```
python tcga_patients.py -p <path/to/yours/netMHCpan4.1> -f <path/to/yours/hg38.fa>
```
