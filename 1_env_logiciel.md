---
title: Préparer l'environnement logiciel
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

# Préparer l'environnement logiciel 🧰

On suppose que conda et mamba sont installés.

## Créer l'environnement conda

Créez un environnement conda :

```bash
conda create -n rnaseq-scere -y
```

Activez ensuite cet environnement :

```bash
conda create -n rnaseq-scere -y
```

## Installer les logiciels nécessaires

Nous allons utiliser [SRA Toolkit](https://github.com/ncbi/sra-tools) pour télécharger les données brutes de séquençage.

L'article indique :

> Raw FASTQ files were aligned to the respective yeast genomes using STAR [78]. Aligned reads were assembled into transcripts, quantified, and normalized using Cufflinks2 [79]. Samples from each yeast time series were normalized together using the CuffNorm feature.

Le fichier *S1 Supporting Information Methods* fournit des précisions supplémentaires :

> Raw FASTQ files from each experiment were aligned to the respective yeast reference genome using STAR [1].

> Reads mapping uniquely to annotated gene features were quantified using HTSeq-count [3].

> Transcript quantification of annotated yeast genes was performed using alignment files output from STAR and Cufflinks2 [4]. Time point samples from the respective yeasts were then normalized together using the CuffNorm feature.

En résumé, nous avons besoin d'installer les outils : `STAR`, `HTSeq-count` et `Cufflinks`. Aucune version de logiciel n'étant spécifiée, nous allons installer la dernière version disponible.

Nous installerons également [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) pour contrôler la qualité des *reads*.

Dans l'environnement conda `rnaseq-scere`, installez tous les logiciels nécessaires :

```bash
mamba install -c conda-forge -c bioconda sra-tools fastqc star htseq cufflinks -y
```

Vérifiez alors les différentes versions des logiciels :

```bash
$ fasterq-dump --version

"fasterq-dump" version 2.11.0
```

```bash
$ fastqc --versio
FastQC v0.11.9
```

```bash
 STAR --version
2.7.10a
```

```bash
$ htseq-count --version
0.13.5
```

```bash
$ cufflinks 2>&1 | head -n 1
cufflinks v2.2.1
```