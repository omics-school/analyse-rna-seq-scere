import os

DATA_DIR = os.environ.get("DATA_DIR", os.getcwd())
BASE_DIR = os.environ.get("BASE_DIR", os.getcwd())

SAMPLES, = glob_wildcards(DATA_DIR + "/reads/{sample}.fastq.gz")
#SAMPLES = ['SRR3405791', 'SRR3405788', 'SRR3405783', 'SRR3405784', 'SRR3405789']

#print(SAMPLES)

onstart:
    print(SAMPLES)

rule make_all:
    input:
        # Quality control
        expand(BASE_DIR + "/reads_qc/{sample}_fastqc.html", sample=SAMPLES),
        # Genome index
        BASE_DIR + "/genome_index/SAindex",
        # BAM sorting and indexing
        expand(BASE_DIR + "/reads_map/{sample}_Aligned.sorted.out.bam", sample=SAMPLES),


rule clean:
    shell:
        "rm -rf reads_qc genome_index reads_map"


rule control_read_quality:
    input:
        DATA_DIR + "/reads/{sample}.fastq.gz"
    output:
        BASE_DIR + "/reads_qc/{sample}_fastqc.html"
    params:
        folder = directory(BASE_DIR + "/reads_qc")
    message:
        "Checking read quality for {input}"
    conda:
        "workflow.yml"
    shell:
        "fastqc {input} "
        "--outdir {params.folder}"


rule index_genome:
    input:
        sequence = DATA_DIR + "/genome/genome.fa",
        annotation = DATA_DIR + "/genome/genes.gtf"
    output:
        index = BASE_DIR + "/genome_index/SAindex"
    params:
        folder = directory(BASE_DIR + "/genome_index")
    message:
        "Indexing reference genome {input.sequence}"
    conda:
        "workflow.yml"
    threads: 
        4
    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.folder} "
        "--genomeFastaFiles {input.sequence} "
        "--sjdbGTFfile {input.annotation} " 
        "--sjdbOverhang 50 "
        "--genomeSAindexNbases 10"


rule map_reads:
    input:
        fastq = DATA_DIR + "/reads/{sample}.fastq.gz",
        index = rules.index_genome.output.index
    output:
        bam = BASE_DIR + "/reads_map/{sample}_Aligned.out.bam"
    params:
        folder = directory(BASE_DIR + "/reads_map"),
        prefix = BASE_DIR + "/reads_map/{sample}_"
    message:
        "Mapping reads {input.fastq} to the reference genome"
    conda:
        "workflow.yml"
    threads: 
        8
    shell:
        "STAR --runThreadN {threads} "
        "--runMode alignReads "
        "--genomeDir {rules.index_genome.params.folder} "
        "--sjdbGTFfile {rules.index_genome.input.annotation} "
        "--readFilesCommand zcat "
        "--readFilesIn {input.fastq} "
        "--outFilterType BySJout "
        "--alignIntronMin 10 "
        "--alignIntronMax 3000 "
        "--outFileNamePrefix {params.prefix} "
        "--outFilterIntronMotifs RemoveNoncanonical "
        "--outSAMtype BAM Unsorted"


rule sort_index_bam:
    input:
        bam = BASE_DIR + "/reads_map/{sample}_Aligned.out.bam"
    output:
        bam = BASE_DIR + "/reads_map/{sample}_Aligned.sorted.out.bam",
        bai = BASE_DIR + "/reads_map/{sample}_Aligned.sorted.out.bam.bai"
    message:
        "Sorting and indexing BAM {input.bam}"
    conda:
        "workflow.yml"
    threads: 
        1
    shell: 
        """
        samtools sort {input.bam} -o {output.bam} 
        samtools index {output.bam}
        """
