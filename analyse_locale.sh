# Le script va s'arrêter :
# - à la première erreur,
# - si une variable n'est pas définie,
# - si une erreur est recontrée dans un pipe.
set -euo pipefail

# Ce script doit être lancé depuis le répertoire de base qui doit exister :
base_dir="/mnt/c/Users/omics/rnaseq_scere"

# Ce script suppose que les échantillons à analyser sont dans le répertoire 'reads'
fastq_dir="${base_dir}/reads"
# et que le génome de référence et ses annotations soit dans le répertoire 'genome'
genome_dir="${base_dir}/genome"

# Références des échantillons à analyser.
# Les numéros d'accession sont entre guillemets et séparés par un espace.
# Faites en sorte que ces numéros correspondent à VOS échantillons dans le répertoire 'reads'.
samples="SRR3405783 SRR3405784 SRR3405788"

# Chemin et nom du fichier contenant le génome de référence.
genome_file="${genome_dir}/genome.fa"
# Chemin et nom du fichier contenant les annotations.
annotation_file="${genome_dir}/genes.gtf"


# On indexe le génome qu'une seule fois.
echo "=============================================================="
echo "Indexer le génome de référence"
echo "=============================================================="
mkdir -p "${base_dir}/genome_index"
STAR --runMode genomeGenerate \
--genomeDir "${base_dir}/genome_index" \
--genomeFastaFiles "${genome_file}" \
--sjdbGTFfile "${annotation_file}" \
--sjdbOverhang 50 \
--genomeSAindexNbases 10


for sample in ${samples}
do
    echo "=============================================================="
    echo "Contrôler la qualité : échantillon ${sample}"
    echo "=============================================================="
    mkdir -p "${base_dir}/reads_qc"
    fastqc "${fastq_dir}/${sample}.fastq.gz" --outdir "${base_dir}/reads_qc"

    echo "=============================================================="
    echo "Aligner les reads sur le génome de référence : échantillon ${sample}"
    echo "=============================================================="
    mkdir -p "${base_dir}/reads_map"
    STAR --runThreadN 1 \
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
    cuffquant --library-type=fr-firststrand "${annotation_file}" \
    "${base_dir}/reads_map/${sample}_Aligned.sortedByCoord.out.bam" \
    --output-dir "${base_dir}/counts/${sample}"
done

echo "=============================================================="
echo "Normaliser les comptages des transcrits"
echo "=============================================================="
cuffnorm --library-type=fr-firststrand "${annotation_file}" \
"${base_dir}/counts"/*/*.cxb --output-dir "${base_dir}/counts"
