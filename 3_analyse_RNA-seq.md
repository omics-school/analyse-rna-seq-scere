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
--sjdbOverhang 50
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

L’indexation du génome n’est à faire qu’une seule fois pour chaque logiciel d’alignement.



