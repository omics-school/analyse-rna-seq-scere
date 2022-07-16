#!/bin/bash

# Script will stop at the first error
set -euo pipefail

# Maximum number of cpus used simultaneously
NJOBS=50

# Load modules
module load snakemake
module load slurm-drmaa

# Paths
export DATA_DIR="/shared/projects/form_2021_29/data/rnaseq_scere"

# Run snakemake with cluster mode
snakemake --profile snakemake_profiles/cluster/ \
    --jobs "${NJOBS}" \
    --drmaa " --mem={resources.mem_mb} --cpus-per-task={threads}"

