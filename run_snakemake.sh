#!/bin/bash

# Script will stop at the first error
set -euo pipefail


# Parameters
NJOBS=50

# Load modules
module load snakemake
module load slurm-drmaa

# --jobs=xx is the maximum number of cpus used simultaneously
#snakemake --cluster-config cluster.yml \
#    --drmaa " --mem={cluster.mem} --cpus-per-task={cluster.cpus}" \
#    --use-conda --jobs "${NJOBS}" --latency-wait 20 --keep-going


snakemake --profile snakemake_profiles/cluster/ \
    --jobs "${NJOBS}" \
    --drmaa " --mem={resources.mem_mb} --cpus-per-task={threads}"

