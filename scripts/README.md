# Stage de fin de master BIBS - Mathilde Bertrand 

(Muséum National d'Histoire Naturelle)
Analyse par ChIP-Seq de l'impact d'une modification ciblée de la chromatine centromérique sur l'ensemble du génome.
Les scripts utilisés pendant le stage pour mettre en place le protocole d'analyse sur les régions uniques et les éléments répétés sont présentés ici. 

## Le dossier Traitement_fastq

Contient les scripts qui ont été utilisés pour le traitement des fichiers fastq et la mise en place du protocole sur les regions uniques et repetees.

```
bash main.sh 0 pour lancer l'analyse qualité
bash main.sh 1 : pour lancer l'analyse sur les régions uniques =>  fait appel à : filtre_strategie1.sh 
bash main.sh 2 : pour lancer l'analyse sur les régions répétées => fait appel à repEnrich.sh 

Remarque : par defaut lanalyse des régions répétées se fait avec RepEnrich2 mais il est possible de la lancer avec RepBase (RepBase.sh)
```

## Le dossier regions_uniques

Ce dossier contient les scripts nécessaires a l'analyse des regions uniques une fois la detection des pics réalisée (cf bash main.sh 1)

### MiseCommuns

Scripts pour la mise en communs des pics entre deux fichiers : distingue les pics chevauchants et les pics uniques

```
python compare_peak.py -f1 file1 -f2 file2 -o output.txt
```

### AnnotGenes

Une fois les pics mis en communs grâce à compare_peak.py, ils peuvent être annootés avec une annotation de gènes : 

```
python annotate_domains_v3.0.py -f MY_doms.txt -t 500 
```

### Genome_decoupe_fenetres

Analyse du genome decoupe en fenetres 

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

Contient les scripts pour generer les fichiers de visualisation sous IGV (IGV.sh) ou avec Circos

## Le dossier Comptages

Les différents scripts qui ont permis de réaliser des comptages

```
Filtre_effects.sh : lors de l'application des filtres pour l'étude des régions uniques, compte le nombre de reads qui restent après application des filtres
gsat_quantifperte.sh : Quantifier la modification de nombre de reads present sur le GSAT du chrs9

```
