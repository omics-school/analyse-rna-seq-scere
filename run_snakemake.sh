#!/bin/bash

# Script will stop at the first error
set -euo pipefail

# Load 
module load snakemake
module load slurm-drmaa

# --jobs=xx is the maximum number of cpus used simultaneously
snakemake --cluster-config cluster.yml --use-conda --drmaa --jobs=50 " --mem={cluster.mem}" 