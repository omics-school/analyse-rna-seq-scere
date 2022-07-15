SAMPLES, = glob_wildcards("reads/{sample}.fastq.gz")


rule make_all:
    input:
        # Quality control
        expand("reads_qc/{sample}_fastqc.html", sample=SAMPLES),
        # Genome index
        "genome_index/SAindex",
        # Reads mapping
        expand("reads_map/{sample}_Aligned.out.bam", sample=SAMPLES),
        # BAM sorting and indexing
        expand("reads_map/{sample}_Aligned.sorted.out.bam", sample=SAMPLES),


rule clean:
    shell:
        "rm -rf reads_qc genome_index reads_map"


rule control_read_quality:
    input:
        "reads/{sample}.fastq.gz"
    output:
        "reads_qc/{sample}_fastqc.html"
    params:
        folder=directory("reads_qc")
    message:
        "Checking read quality for {input}"
    conda:
        "envs/workflow.yml"
    shell:
        "fastqc {input} "
        "--outdir {params.folder}"


rule index_genome:
    input:
        sequence="genome/genome.fa",
        annotation="genome/genes.gtf"
    output:
        index="genome_index/SAindex"
    params:
        folder=directory("genome_index")
    message:
        "Indexing reference genome {input.sequence}"
    conda:
        "envs/workflow.yml"
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
        fastq="reads/{sample}.fastq.gz",
        index=rules.index_genome.output.index
    output:
        bam="reads_map/{sample}_Aligned.out.bam"
    params:
        folder=directory("reads_map"),
        prefix="reads_map/{sample}_"
    message:
        "Mapping reads {input.fastq} to the reference genome"
    conda:
        "envs/workflow.yml"
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
        bam="reads_map/{sample}_Aligned.out.bam"
    output:
        bam="reads_map/{sample}_Aligned.sorted.out.bam",
        bai="reads_map/{sample}_Aligned.sorted.out.bam.bai"
    message:
        "Sorting and indexing BAM {input.bam}"
    conda:
        "envs/workflow.yml"
    threads: 
        1
    shell: 
        """
        samtools sort {input.bam} -o {output.bam} 
        samtools index {input.bam}
        """
