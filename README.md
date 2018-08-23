# Stage de fin de master BIBS - Mathilde Bertrand 

(Muséum National d'Histoire Naturelle)
Analyse par ChIP-Seq de l'impact d'une modification ciblée de la chromatine centromérique sur l'ensemble du génome.
Les scripts utilisés pendant le stage pour mettre en place le protocole d'analyse sur les régions uniques et les éléments répétés et pour générer les résultats présents dans le rapport de stage sont présentés ici. 

## Le dossier Traitement_fastq

Contient les scripts qui ont été utilisés pour le traitement des fichiers fastq et la mise en place du protocole sur les regions uniques et repetees.

```
bash main.sh 0 pour lancer l'analyse qualité
bash main.sh 1 : pour lancer l'analyse sur les régions uniques =>  fait appel à : filtre_strategie1.sh 
bash main.sh 2 : pour lancer l'analyse sur les régions répétées => fait appel à repEnrich.sh 

Remarques : 
- par defaut lanalyse des régions répétées se fait avec RepEnrich2 mais il est possible de la lancer avec RepBase (RepBase.sh)
- lanalyse R des regions repetees sur R est directement lancee, ce qui nest pas le cas pour les regions uniques
```

## Le dossier regions_uniques

Ce dossier contient les scripts nécessaires a l'analyse des regions uniques une fois la detection des pics réalisée (cf bash main.sh 1) : 

### MiseCommuns

Scripts pour la mise en communs des pics entre deux fichiers : distingue les pics chevauchants et les pics uniques

```
python compare_peak.py -f1 file1 -f2 file2 -o output.txt
f1 et f2 sont les fichiers de pics
```

### AnnotGenes

Une fois les pics mis en communs grâce à compare_peak.py, ils peuvent être annotés avec une annotation de gènes au format gff : 

```
python annotate_domains_v3.0.py -f <INPUT> -a <annot> -e <extension> -t <THRESHOLD> -o <OUTPUT> 
```

### Genome_decoupe_fenetres

Analyse du genome decoupe en fenetres de 55, 25, 15 et 5 kb : decoupage du genome en fenêtres, comptage du nombre de reads dans chaque fenetres puis ACP pour determiner des liens entre les échantillons.

 ```
 Pour lancer l'analyse : 
 bash decoupage.sh
 ```

### R
Scripts R pour l'analyse sur les régions uniques

```
Violin.R : script permettant de réaliser les violin plot pour comparer la taille des pics
Manhantann.R : script permettant de réaliser les Manhantann plot pour choisir le seuil de FDR
Quantif-modifIntensite.sh : Quantifier la modification de nombre de reads present sur le GSAT du chrs9
```


## Le dossier regions_repetees

Contient les scripts d'analyses statistique R de Repenrich2 (repEnrich.R) et de RepBase (RepBase.R)



## Le dossier Visualisation

Contient les scripts pour générer les fichiers de visualisation sous IGV (IGV.sh) ou avec Circos

## Le dossier Comptages

Les différents scripts qui ont permis de réaliser des comptages

```
Filtre_effects.sh : lors de l'application des filtres pour l'étude des régions uniques, compte le nombre de reads qui restent après application des filtres
gsat_quantifperte.sh : Quantifier la modification de nombre de reads present sur le GSAT du chrs9

```

## Architecture des dossiers  : 

Obligatoire : la présence des dossiers : raw_data/labo_curie/fastq avec les différents fichiers au format fastq
genome_ref : contenant le génome de référence mm10 et le génome rodent de RepBase au format fasta
Tools : contenant les outils fastq et bigwigtobedgraphe, RepEnrich2, IGV, circos
analysis
Il faut également que les outils suivants soient installés : bowtie2, samtools, bedtools, deeptools, R (avec les packages FactoMineR, DeSeq2), epic, (macs2)

Les dossiers Trimme (les fichiers fastq trimés des 5 premières paires de bases), FASTQC (contenant les éléments d'analyse qualité), mapping (contenant les différents fichiers d'alignments et les fichiers bam filtrés), et PICS (contenant les résultats de la détection de pics), Repenrich, repBaseMapping et Fenetres (analyse du genome decoupes en fenetres) sont crées au moment de l'analyse.

![alt text](https://github.com/MathildeBertrand/STAGE_M2/blob/master/Organisation_dossiers.png)