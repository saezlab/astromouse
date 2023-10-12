#!env bash

echo $CONDA_PREFIX

mkdir -p logs

(test -f logs/celloracle.post-deploy.log || rm logs/celloracle.post-deploy.log)

$CONDA_PREFIX/bin/Rscript workflow/envs/celloracle.R >> logs/celloracle.post-deploy.log 2>&1
