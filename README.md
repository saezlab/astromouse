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

### Run pipeline
Snakemake takes care of installing all dependencies to the environments. Set them up by running:
```
snakemake --conda-create-envs-only
```

You can see which steps will be executed by running:
```
snakemake -n
```

Run all steps necessary to reproduce the plots:
```
snakemake --use-conda -c 4 #will run on 4 local cores
```

For a cluster with SLURM as job manager:
```
snakemake --profile config/slurm -j 1 #will submit only one job at a time
```

More information about how to use and execute snakemake can be found in the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)



