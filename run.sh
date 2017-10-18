#!/usr/bin/env bash
echoerr() { printf "%s\n" "$*" >&2; }

# Load the conda environment
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echoerr "Load the conda environment ..."
    source /opt/conda/bin/activate /gscmnt/gc7201/dinglab/medseq/lwang/conda_root/envs/gatk_wxs_cnv
fi

# Run snakemake via bsub
echoerr "Launch the LSF/bsub master job to run snakemake ..."
bsub -a "docker(lbwang/dailybox)" -q research-hpc -N -o snakemake_job.log \
    env SHELL="/bin/bash" snakemake --timestamp --rerun-incomplete --nolock \
    --jobs 100 \
    --cluster "./bsub_submitter.py {dependencies} lsf_logs" \
    --cluster-config bsub_config.json \
    -p all
