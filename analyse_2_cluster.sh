#!/bin/bash

#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-3            # limit to 4 samples. Use --array=0-49 for all 50 samples.

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# chargement des modules nécessaires
module load star/2.7.9a


# répertoire de base (le répertoire depuis lequel vous lancez le script)
base_dir="$PWD"
# répertoire contenant les données
data_dir="/shared/projects/form_2021_29/data/rnaseq_scere"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="${data_dir}/genome"
# chemin et nom du fichier contenant le génome de référence.
genome_file="${genome_dir}/genome.fa"
# chemin et nom du fichier contenant les annotations
annotation_file="${genome_dir}/genes.gtf"
# répertoire contenant les fichiers .fastq.gz
fastq_dir="${data_dir}/reads"
# liste de tous les fichiers .fastq.gz
fastq_files=(${fastq_dir}/*fastq.gz)
# extraction de l'identifiant de l'échantillon
# à partir du nom de fichier : /shared/projects/form_2021_29/data/rnaseq_tauri/reads/SRR3405783.fastq.gz
# on extrait : SRR3405783
sample=$(basename -s .fastq.gz "${fastq_files[$SLURM_ARRAY_TASK_ID]}")


echo "=============================================================="
echo "Contrôler la qualité : échantillon ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/reads_qc"
fastqc "${fastq_dir}/${sample}.fastq.gz" --outdir "${base_dir}/reads_qc"


echo "=============================================================="
echo "Aligner les reads sur le génome de référence : échantillon ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/reads_map"
STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
--runMode alignReads \
--genomeDir "${base_dir}/genome_index" \
--sjdbGTFfile "${annotation_file}" \
--readFilesCommand zcat \
--readFilesIn "${fastq_dir}/${sample}.fastq.gz" \
--outFilterType BySJout \
--alignIntronMin 10 \
--alignIntronMax 3000 \
--outFileNamePrefix "${base_dir}/reads_map/${sample}_" \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM SortedByCoordinate

echo "=============================================================="
echo "Indexer les reads alignés : échantillon ${sample}"
echo "=============================================================="
samtools index "${base_dir}/reads_map/${sample}_Aligned.sortedByCoord.out.bam"

echo "=============================================================="
echo "Compter les reads : échantillon ${sample}"
echo "=============================================================="
mkdir -p "${base_dir}/counts/${sample}"
htseq-count --order=pos --stranded=reverse \
--mode=intersection-nonempty \
"${base_dir}/reads_map/${sample}_Aligned.sortedByCoord.out.bam" \
"${annotation_file}" > "${base_dir}/counts/${sample}/count.txt"

echo "=============================================================="
echo "Compter les transcrits : échantillon ${sample}"
echo "=============================================================="
cuffquant --num-threads "${SLURM_CPUS_PER_TASK}" \
--library-type=fr-firststrand "${annotation_file}" \
"${base_dir}/reads_map/${sample}_Aligned.sortedByCoord.out.bam" \
--output-dir "${base_dir}/counts/${sample}"