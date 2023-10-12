#!env bash
conda activate -p $CONDA_PREFIX

mkdir -p logs

(test -f logs/misty.post-deploy.log || rm logs/misty.post-deploy.log)

$CONDA_PREFIX/bin/Rscript workflow/envs/misty.R >> logs/misty.post-deploy.log 2>&1
