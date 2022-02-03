---
title: Analyser les données RNA-ses
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

# Analyser les données RNA-seq 💻

## Préparer l'environnement

Si cela n'est pas déjà fait, activez l’environnement conda *rnaseq-scere* :

```bash
conda activate rnaseq-scere
```

Déplacez-vous ensuite dans le répertoire contenant les répertoires `reads` et `genome` de *S. cerevisiae*. Vous devriez obtenir l'arborescence suivante :

```bash
$ tree
.
├── genome
│   ├── genes.gtf
│   └── genome.fa
└── reads
    ├── SRR3405783.fastq.gz
    ├── SRR3405784.fastq.gz
    ├── SRR3405788.fastq.gz
    ├── SRR3405789.fastq.gz
    └── SRR3405791.fastq.gz
```

## Analyser manuellement un échantillon

Choississez un échantillon parmi ceux téléchargés dans le répertoire `reads` :

```bash
$ tree reads
reads
├── SRR3405783.fastq.gz
├── SRR3405784.fastq.gz
├── SRR3405788.fastq.gz
├── SRR3405789.fastq.gz
└── SRR3405791.fastq.gz
```

Par exemple : `SRR3405783.fastq.gz`


### Contrôler la qualité des reads

Créez le répertoire `reads_qc` qui va contenir les fichiers produits par le contrôle qualité des fichiers *fastq.gz* :

```bash
mkdir -p reads_qc
```

Lancez FastQC avec la commande :

```bash
fastqc reads/SRR3405783.fastq.gz --outdir reads_qc
```

FastQC va produire deux fichiers (un fichier avec l’extension `.html` et un autre avec l’extension `.zip`) dans le répertoire `reads_qc`. Si par exemple, vous avez analysé le fichier `reads/SRR3405783.fastq.gz`, vous obtiendrez les fichiers `reads_qc/SRR3405783_fastqc.html` et `reads_qc/SRR3405783_fastqc.zip`.

Ouvrez le fichier `.html` ainsi créé avec Firefox. Analysez le rapport créé par FastQC.


### Indexer le génome de référence

L’indexation du génome de référence est une étape indispensable pour accélérer l’alignement des reads sur le génome. Elle consiste à créer un annuaire du génome de référence.

Créez le répertoire `genome_index` qui contiendra les index du génome de référence :

```bash
mkdir -p genome_index
```

Lancez l’indexation du génome de référence :

```bash
STAR --runMode genomeGenerate \
--genomeDir genome_index \
--genomeFastaFiles genome/genome.fa \
--sjdbGTFfile genome/genes.gtf \
--sjdbOverhang 50 \
--genomeSAindexNbases 10
```

L'aide de STAR pour l'option `--sjdbOverhang` indique :

```
sjdbOverhang                            100
    int>0: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)
```

Cela signifie que cette option doit être égale à la longueur maximale des *reads* - 1. Le fichier *S1 Supporting Information Methods* mentionne :

> S. cerevisiae total mRNA was prepared in libraries of stranded 50 base-pair single-end reads and multiplexed at 10 time point samples per sequencing lane.

Les *reads* obtenus devraient donc a priori être constitué de 50 bases. On peut le vérifier en affichant les premières lignes d'un fichier *.fastq.gz*, par exemple `reads/SRR3405783.fastq.gz` :

```bash
$ zcat reads/SRR3405783.fastq.gz | head
@SRR3405783.1 3NH4HQ1:254:C5A48ACXX:1:1101:1135:2105/1
GGTTGAANGGCGTCGCGTCGTAACCCAGCTTGGTAAGTTGGATTAAGCACT
+
?8?D;DD#2<C?CFE6CGGIFFFIE@DFF<FFB===C7=F37@C)=DE>EA
@SRR3405783.2 3NH4HQ1:254:C5A48ACXX:1:1101:1210:2111/1
GTTTCTGTACTTACCCTTACCGGCTCTCAATTTCTTGGACTTCAAGACCTT
+
;?@DDD>BFF?F>GE:CFFFFBFFIFFFIIIEIFFI9BFDCFFC<FFF>FF
@SRR3405783.3 3NH4HQ1:254:C5A48ACXX:1:1101:1408:2208/1
CTGCAGACAAGGCATCTCCTCTCAAGGCCAAATGACGTTGGTCCAACAGTA
```

En réalité, les *reads* ont une longueur de 51 bases et non pas 50 comme supposé.

Le paramètre `--sjdbOverhang` vaut donc 50 (51 - 1).

*Remarque 1 : La commande `zcat` est particulière car elle affiche le contenu d'un fichier compressé (ici un fichier .gz) en le décompressant à la volée. La commande `cat` (sans le `z`) n'aurait pas permis une telle manipulation.*

Vérifiez également que la longueur des *reads* est bien 51 en consultant les résultats du contrôle qualité fournis par FastQC.

Enfin, le paramètre `--genomeSAindexNbases 10` est conseillé par STAR si on le lance sans :

> !!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=12157105, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 10

Nous vous rappelons que l’indexation du génome n’est à faire qu’une seule fois pour chaque logiciel d’alignement.


### Aligner les *reads* sur le génome de référence

Le fichier *S1 Supporting Information Methods* précise la commande utilisée pour l'alignement :

```bash
STAR --runThreadN 1 --runMode alignReads --genomeDir
path_to_yeast_genome_build --sjdbGTFfile path_to_yeast_transcriptome_gtf
--readFilesIn sample.fastq --outFilterType BySJout --alignIntronMin 10 --
alignIntronMax 3000 --outFileNamePrefix ./STAR_out/ --
outFilterIntronMotifs RemoveNoncanonical
```

Les alignements des *reads* seront stockés dans le répertoire `reads_map` :

```bash
mkdir -p reads_map
```

Avec nos chemins de fichiers et quelques adaptations, la commande d'alignement devient :

```bash
STAR --runThreadN 1 \
--runMode alignReads \
--genomeDir genome_index \
--sjdbGTFfile genome/genes.gtf \
--readFilesCommand zcat \
--readFilesIn reads/SRR3405783.fastq.gz \
--outFilterType BySJout \
--alignIntronMin 10 \
--alignIntronMax 3000 \
--outFileNamePrefix reads_map/SRR3405783_ \
--outFilterIntronMotifs RemoveNoncanonical \
--outSAMtype BAM Unsorted
```

**Remarques** :

- L'article orginal de STAR « [STAR: ultrafast universal RNA-seq aligner](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) » (*Bioinformatics*, 2013) précise que :
    > STAR’s default parameters are optimized for mammalian genomes. Other species may require significant modifications of some alignment parameters; in particular, the maximum and minimum intron sizes have to be reduced for organisms with smaller introns.
    
    Ceci explique pourquoi les options `--alignIntronMin 10` et `--alignIntronMax 3000` ont été adaptées pour le génome de la levure *S. cerevisiae*.

- L'option `--readFilesCommand zcat` n'était pas présente dans la commande fournie en *Supporting information*. Nous l'avons ajoutée car les fichiers contenant les *reads* (*.fastq.gz*) sont compressés et il faut demander explicitement à STAR de le prendre en charge. Pensez à consultez toujour la [documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) de l'outil que vous utilisez !

Lancez l'alignement avec STAR et vérifiez que tout se déroule sans problème.

### Compter les *reads* et les transcrits

Le fichier *S1 Supporting Information Methods* précise les commandes utilisées pour le comptage des *reads* :

```bash
htseq-count --order=pos --stranded=reverse --mode=intersection-nonempty
sample.aligned.sorted.sam path_to_yeast_transcriptome_gtf > sample.txt
```

et celui des transcrits :

```bash
cuffquant --library-type=fr-firststrand path_to_yeast_transcriptome_gtf
sample.aligned.sorted.sam
```

Enfin, la normalisation des comptages des transcrits :

```bash
cuffnorm --library-type=fr-firststrand path_to_yeast_transcriptome_gtf
*.cxb
```

Nous stockons les fichiers de comptage dans le répertoire `counts/SRR3405783` :

```bash
mkdir -p counts/SRR3405783
```

Les étapes de tri et d'indexation des *reads* alignés ne sont pas explicitement mentionnées dans les *Supporting Informations* mais elles sont cependant nécessaires pour `HTSeq` :

```bash
samtools sort reads_map/SRR3405783_Aligned.out.bam -o reads_map/SRR3405783_Aligned.sorted.out.bam
samtools index reads_map/SRR3405783_Aligned.sorted.out.bam
```

Remarque : le tri des *reads* peut, a priori, se faire directement avec STAR en utilisant l'option `--outSAMtype BAM SortedByCoordinate`. Cependant, les tests réalisés sur le cluster ont montré qu'il y avait parfois des soucis avec cette option. Nous préférons donc utiliser `samtools` pour le tri des *reads*.


La commande pour compter les *reads* devient alors :

```bash
htseq-count --order=pos --stranded=reverse \
--mode=intersection-nonempty \
reads_map/SRR3405783_Aligned.sorted.out.bam genome/genes.gtf > counts/SRR3405783/count_SRR3405783.txt
```

Puis celle pour compter les transcrits :

```bash
cuffquant --library-type=fr-firststrand genome/genes.gtf \
reads_map/SRR3405783_Aligned.sorted.out.bam \
--output-dir counts/SRR3405783
```

Remarque :

- Par défaut, `cuffquant` écrit un fichier `abundances.cxb`.
- Nous ajoutons l'option `--output-dir counts/SRR3405783` pour indiquer où stocker les résultats produits par `cuffquant` (voir [documentation](http://cole-trapnell-lab.github.io/cufflinks/cuffquant/)). Cela nous permet de distinguer les résultats obtenus à partir de différents fichiers *.fastq.gz*.


Enfin, on normalise les comptages des transcrits :

```bash
cuffnorm --library-type=fr-firststrand genome/genes.gtf \
counts/*/*.cxb --output-dir counts
```

Remarques : 

- Dans le cas présent, cette normalisation va échouer car nous n'avons aligné et quantifié qu'un seul fichier *.fastq.gz*. Cette étape sera par contre pertinente lorsque plusieurs fichiers *.fastq.gz* seront traités.
- L'option `--output-dir counts` indique où stocker les fichiers produits par `cuffnorm` (voir [documentation](http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/)).


## Analyser automatiquement 3 échantillons

Vérifiez que vous êtes bien dans le répertoire `/mnt/c/Users/omics/rnaseq_scere`. Assurez-vous également que vous avez préparé les données correctement, notamment les répertoires `reads` et `genome` :

```bash
$ tree
.
├── genome
│   ├── genes.gtf
│   └── genome.fa
└── reads
    ├── SRR3405783.fastq.gz
    ├── SRR3405784.fastq.gz
    ├── SRR3405788.fastq.gz
    ├── SRR3405789.fastq.gz
    └── SRR3405791.fastq.gz
```

Téléchargez le script `analyse_locale.sh` qui analyse 3 échantillons :

```bash
wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq-scere/master/analyse_locale.sh
```

Lancez ensuite le script d'analyse :

```bash
bash analyse_locale.sh
```

L'analyse devrait prendre plusieurs dizaines de minutes.

Vérifiez régulièrement votre terminal qu'aucune erreur n'apparaît.

Le fichier qui contient le comptage normalisé des transcrits est `counts/genes.count_table`.

## Analyser automatiquement 50 échantillons sur un cluster

Connectez-vous sur le cluster de l'IFB.

Déplacez-vous dans votre répertoire de travail :

```bash
cd "/shared/projects/form_2021_29/$USER"
```

puis créez le répertoire `rnaseq_scere` :

```bash
mkdir -p rnaseq_scere
```

et déplacez-vous à l'intérieur :

```bash
cd "/shared/projects/form_2021_29/$USER/rnaseq_scere"
```

Vérifiez que toutes les données sont bien dans `/shared/projects/form_2021_29/data/rnaseq_scere` :

```bash
tree /shared/projects/form_2021_29/data/rnaseq_scere
```

L'analyse complète des données RNA-seq de *Saccharomyces cerevisiae* (50 fichiers *.fastq.gz*) va se faire en 3 étapes :

1. Indexer le génome. Cette étape est à réaliser une seule fois.
2. Aligner les *reads* sur le génome de référence. Cette étape est particulière car elle sera distribuée sur autant de jobs qu'il y a de fichiers *.fastq.gz* à analyser.
3. Normaliser les transcrits. Cette étape est réalisée en dernier et nécessite que tous les jobs de l'étape 2 soient terminés.

### Indexer le génome de référence

Téléchargez le script `analyse_1_cluster.sh` :

```bash
wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq-scere/master/analyse_1_cluster.sh
```

puis lancez ensuite ce script :

```bash
sbatch -A form_2021_29 analyse_1_cluster.sh
```

Vérifiez que votre job est bien lancé avec la commande :

```bash
squeue -u $USER
```

Le fichier `slurm-JOBID.out` est également créé et contient les sorties du script. Pour consulter son contenu, tapez :

```bash
cat slurm-JOBID.out
```

avec `JOBID` le numéro de votre job.

Suivez également en temps réel l'exécution de votre job avec la commande :

```bash
watch sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j JOBID
```

avec `JOBID` le numéro de votre job.


## Aligner les *reads* sur le génome de référence

Téléchargez le script `analyse_2_cluster.sh` :

```bash
wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq-scere/master/analyse_2_cluster.sh
```

puis lancez ensuite ce script :

```bash
sbatch -A form_2021_29 analyse_2_cluster.sh
```

Suivez l'évolution de votre job avec les commandes habituelles.


## Normaliser les comptages des transcrits

Téléchargez le script `analyse_3_cluster.sh` :

```bash
wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq-scere/master/analyse_3_cluster.sh
```

puis lancez ensuite ce script :

```bash
sbatch -A form_2021_29 analyse_3_cluster.sh
```

Suivez l'évolution de votre job avec les commandes habituelles.

Les données de comptage agrégées et normalisées se trouvent dans le fichier :

```bash
/shared/projects/form_2021_29/$USER/rnaseq_scere/counts/counts/genes.count_table
```

## Bonus : aggréger les données produites par HTSeq-count

Si vous analysez plusieurs échantillons, vous souhaiterez peut-être aggréger tous les fichiers produits par HTSeq-count.

Les instructions Bash suivantes pourront vous y aider :

```bash
for name in counts/*/*.txt
do
    echo "${name}"
    # On récupère le nom de l'échantillon
    sample="$(basename -s .txt ${name} | sed 's/count_//g' )"
    # On stocke dans un fichier temporaire pour chaque échantillon
    # un entête avec "gene" et le nom de l'échantillon
    echo -e "${sample}" > "count_${sample}_tmp.txt"
    # On copie le contenu du fichier de comptage dans ce fichier temporaire
    cut -f2 "${name}" >> "count_${sample}_tmp.txt"
done

# On récupère les noms des gènes
echo "gene" > genes.txt
cut -f1 "${name}" >> genes.txt

# On fusionne tous les fichiers
paste genes.txt *tmp.txt > count_all.txt

# On supprime les lignes qui débutent par '__'
# et qui ne sont pas utiles
grep -v "^__" count_all.txt > count_all_clean.txt

# On supprime les fichiers temporaires
rm -f genes.txt *tmp.txt count_all.txt
```

