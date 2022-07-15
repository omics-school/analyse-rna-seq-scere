#!/bin/bash

# Script will stop at the first error
set -euo pipefail

# Load 
module load snakemake
module load slurm-drmaa

snakemake --cluster-config cluster.yml --use-conda --drmaa --jobs=20 " --mem={cluster.mem}" 