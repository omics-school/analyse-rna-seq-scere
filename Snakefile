import os

DATA_DIR = os.environ.get("DATA_DIR", os.getcwd())

SAMPLES = ["SRR3405791", "SRR3405788", "SRR3405783", "SRR3405784", "SRR3405789"]
#SAMPLES, = glob_wildcards(DATA_DIR + "/reads/{sample}.fastq.gz")

workdir: os.environ.get("WORK_DIR", os.getcwd())


onstart:
    print(f"Working directory:\n{os.getcwd()}")
    print(f"Samples:\n{SAMPLES}")


rule make_all:
    input:
        # Quality control
        expand("reads_qc/{sample}_fastqc.html", sample=SAMPLES),
        # Genome index
        "genome_index/SAindex",
        # BAM sorting and indexing
        expand("reads_map/{sample}_Aligned.sorted.out.bam", sample=SAMPLES),


rule control_read_quality:
    input:
        DATA_DIR + "/reads/{sample}.fastq.gz"
    output:
        "reads_qc/{sample}_fastqc.html"
    params:
        folder = directory("reads_qc")
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
        index =  "genome_index/SAindex"
    params:
        folder = directory("genome_index")
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
        bam = "reads_map/{sample}_Aligned.out.bam"
    params:
        folder = directory("reads_map"),
        prefix = "reads_map/{sample}_"
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
        bam = "reads_map/{sample}_Aligned.out.bam"
    output:
        bam = "reads_map/{sample}_Aligned.sorted.out.bam",
        bai = "reads_map/{sample}_Aligned.sorted.out.bam.bai"
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
