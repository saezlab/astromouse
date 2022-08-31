# Astromouse
Multiome and visium data from the brain and heart of mice that went to space.

This pipeline uses `Snakemake` to ensure reproducibility.

## Installation
Clone repo:
```
git clone git@github.com:saezlab/astromouse.git
cd astromouse
```

Install `mamba` (this might take a while) to install packages faster:
```
conda install -n base -c conda-forge mamba
```

Then create a new enviroment specific for `Snakemake`:
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

