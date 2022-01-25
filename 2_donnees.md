---
title: Préparer les données
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

# Préparer les données 🗃️

L'article orginale publié en 2016 par [Kelliher *et al.*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006453) indique dans la rubrique *Data Availability* :

> RNA-Sequencing gene expression data from this manuscript have been submitted to the NCBI Gene Expression Omnibus (GEO; http://www.ncbi.nlm.nih.gov/geo/) under accession number GSE80474.

Le numéro du projet qui nous intéresse est donc : **GSE80474**

## Données de séquençage

### Méthode 1 : SRA Run Selector

Dans l'outil [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/), entrez l'identifiant du projet : GSE80474. 

Un total de 74 *runs* sont disponibles. Cliquez alors sur le bouton gris *Metadata* correspondant au total.

Téléchargez le fichier `SraRunTable.txt` proposé. Il s'agit d'un fichier CSV, c'est-à-dire d'un fichier tabulé avec des colonnes séparées par des virgules. Ouvrez-le ensuite avec Microsoft Excel ou LibreOffice Calc.

Sélectionnez tous les échantillons correspondants à *Saccharomyces cerevisiae* (les 50 premiers) et copiez / collez la colonne *Run* dans un fichier texte (`runs_scere.txt`)

On peut également réaliser cette sélection en ligne de commande :

```bash
grep "Saccharomyces cerevisiae" SraRunTable.txt | cut -d"," -f1 > runs_scere.txt
```

pour avoir les 5 premiers échantillons uniquement :

```bash
grep "Saccharomyces cerevisiae" SraRunTable.txt | cut -d"," -f1 | head -n 5 > runs_scere_small.txt
```

Téléchargez les fichiers fastq :

```bash
mkdir -p reads
for sample in $(cat runs_scere_small.txt)
do 
    echo ${sample}
    fasterq-dump --progress --outdir reads ${sample}
done
```

```bash
$ du -csh reads/*
4,1G    reads/SRR3405783.fastq
4,7G    reads/SRR3405784.fastq
4,3G    reads/SRR3405785.fastq
4,0G    reads/SRR3405786.fastq
4,0G    reads/SRR3405787.fastq
21G     total
```

Compressez les fichiers fastq :

```bash
gzip reads/*
```

```bash
$ du -csh reads/*
864M    reads/SRR3405783.fastq.gz
983M    reads/SRR3405784.fastq.gz
910M    reads/SRR3405785.fastq.gz
834M    reads/SRR3405786.fastq.gz
835M    reads/SRR3405787.fastq.gz
4,4G    total
```


### Méthode 2 : SRA Explorer

Le numéro du projet GSE80474 commence par les lettres `GSE` ce qui nous indique que c'est un projet initialement déposé dans la base de données GEO. Cette base n'étant pas toujours bien prise en charge par l'outil sra-explorer, nous allons tout d'abord récupérer sur le site SRA Run Selector l'identifiant *BioProject* correspondant.

Sur le site [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/) :

1. Entrez l'identifiant du projet : GSE80474.
2. Cliquez sur le bouton bleu *Search*.
3. Dans la rubrique *Common Fields*, récupérez l'identifiant *BioProject* : PRJNA319029

C'est avec cet identifiant BioProject que nous allons récupérer les données.

Sur le site [SRA EXplorer](https://sra-explorer.info/) :

1. Indiquez le numéro du *BioProject*, ici PRJNA319029, puis cliquez sur le petite loupe pour lancer la recherche.
1. Vous obtenez ensuite 74 réponses qui correspondent aux différents fichiers / échantillons.
1. Affinez les réponses en tapant « Scerevisiae » dans le champ « Filter results: ». vous devriez obtenir 50 résultats.
1. Sélectionnez tous les résultats en cliquant sur le case voide à droite de *Title*.
1. Cliquez sur le bouton « Add 50 to collection ».
1. Cliquez ensuite en haut à droite sur le bouton « 50 saved datasets ».
1. Cliquez enfin sur « Bash script for downloading FastQ files ».
1. Téléchargez le script qui vous permettra de télécharger tous les fichiers fastq (`sra_explorer_fastq_download.sh`).

Voici les 5 premières lignes du script téléchargé :

```bash
$ head -n 5 sra_explorer_fastq_download.sh
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/009/SRR3405789/SRR3405789.fastq.gz -o SRR3405789_GSM2128026_Scerevisiae_YEPD_aF_30min_Saccharomyces_cerevisiae_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/001/SRR3405791/SRR3405791.fastq.gz -o SRR3405791_GSM2128028_Scerevisiae_YEPD_aF_40min_Saccharomyces_cerevisiae_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/004/SRR3405784/SRR3405784.fastq.gz -o SRR3405784_GSM2128021_Scerevisiae_YEPD_aF_5min_Saccharomyces_cerevisiae_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/003/SRR3405783/SRR3405783.fastq.gz -o SRR3405783_GSM2128020_Scerevisiae_YEPD_aF_0min_Saccharomyces_cerevisiae_RNA-Seq.fastq.gz
```

C'est bien un script Bash car la première ligne est `#!/usr/bin/env bash`. Ensuite, chaque ligne qui débute par `curl` télécharge un fichier .fastq.gz. La syntaxe de la commande `curl` est la suivante :

```bash
curl -L ADRESSE-DU-FICHIER-À-TÉLÉCHARGER -o NOM-DU-FICHIER-SUR-LE-DISQUE-LOCAL
```

Nous aimerions modifier ce script pour faire en sorte que 

1. Le nom du fichier enregistré localement ne contienne que le numéro d'accession du fichier, tel que présent sur les serveurs de SRA (par exemple : `SRR3405789`) et pas les métadonnnées associées (par exemple : `_GSM2128026_Scerevisiae_YEPD_aF_30min_Saccharomyces_cerevisiae_RNA-Seq.fastq.gz`). Il faut remplacer l'option `-o` par `-O` (sans argument).
2. Tous les fichiers soient enregistrés dans le même répertoire (par exemple `reads`). Il faut alors ajouter l'option `--output-dir` avec l'arguement `reads`.

Nous utilisons ici la commande `sed` qui peut modifier à la volée les lignes d'un fichier :

```bash
sed -E 's/-o .*/-O --output-dir reads/' sra_explorer_fastq_download.sh  > sra_explorer_fastq_download_2.sh
```

Voici les 5 premières lignes du script `sra_explorer_fastq_download_2.sh` : 

```bash
$ head -n 5 sra_explorer_fastq_download_2.sh
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/009/SRR3405789/SRR3405789.fastq.gz -O --output-dir reads
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/001/SRR3405791/SRR3405791.fastq.gz -O --output-dir reads
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/004/SRR3405784/SRR3405784.fastq.gz -O --output-dir reads
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/003/SRR3405783/SRR3405783.fastq.gz -O --output-dir reads
```

Pour vérifier que notre script fonctionne, nous allons télécharger les 5 premiers fichiers *.fastq.gz*. Pour cela, créez un script intermédiaire en ne sélectionnant que les 6 premières lignes du script de téléchargement :

```bash
head -n 6 sra_explorer_fastq_download_2.sh > sra_explorer_fastq_download_2_small.sh
```

Puis téléchargez les données :

```bash
mkdir -p reads
bash sra_explorer_fastq_download_2_small.sh
```

Patientez quelques minutes que le téléchargement se termine, puis calculez l'espace occupé par les données :

```bash
$ du -csh reads/*
776M    reads/SRR3405783.fastq.gz
883M    reads/SRR3405784.fastq.gz
831M    reads/SRR3405788.fastq.gz
899M    reads/SRR3405789.fastq.gz
901M    reads/SRR3405791.fastq.gz
4,2G    total
```

## Génome de référence et annotations

On trouve dans le fichier *S1 Supporting Information Methods* la desciption du génome de *S. cerevisiae* utilisé :

> The S. cerevisiae S288C genome (Ensembl build R64-1-1) was downloaded from Illumina iGenomes on March 2, 2016 (https://support.illumina.com/sequencing/sequencing_software/igenome.html).

Téléchargez le fichier concerné (*Ensembl build R64-1-1*) et décompressez l'archive :

```bash
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz
tar zxvf Saccharomyces_cerevisiae_Ensembl_R64-1-1.tar.gz
```

Récupérez ensuite les fichiers contenant le génome et les annotations :
```bash
mkdir -p genome
cp Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.gtf genome
cp Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/WholeGenomeFasta/genome.fa genome
```

Supprimez enfin le répertoire `Saccharomyces_cerevisiae` qui ne nous intéresse plus :

```bash
rm -rf Saccharomyces_cerevisiae
```

