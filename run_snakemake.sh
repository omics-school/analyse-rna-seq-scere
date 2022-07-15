#!/bin/bash

# Script will stop at the first error
set -euo pipefail


# Parameters
NJOBS=30

# Load modules
module load snakemake
module load slurm-drmaa

# --jobs=xx is the maximum number of cpus used simultaneously
snakemake --cluster-config cluster.yml \
    --drmaa " --mem={cluster.mem}" \
    --use-conda \
    --jobs="${NJOBS}"