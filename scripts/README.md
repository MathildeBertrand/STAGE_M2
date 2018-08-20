# Stage de fin de master BIBS - Mathilde Bertrand 

(Muséum National d'Histoire Naturelle)
Analyse par ChIP-Seq de l'impact d'une modification ciblée de la chromatine centromérique sur l'ensemble du génome.
Les scripts utilisés pendant le stage pour mettre en place le protocole d'analyse sur les régions uniques et les éléments répétés sont présentés ici. 

##Le dossier Traitement_fastq

Contient les scripts qui ont été utilisés pour le traitement des fichiers fastq et la mise en place du protocole sur les regions uniques et repetees.

```
bash main.sh 0 pour lancer l'analyse qualité
bash main.sh 1 : pour lancer l'analyse sur les régions uniques
bash main.sh 2 : pour lancer l'analyse sur les régions répétées 
``

##Le dossier regions_uniques

Ce dossier contient les scripts nécessaires a l'analyse des regions uniques une fois la detection des pics réalisée (cf bash main.sh 1)

##Le dossier regions_repetees

Contient les scripts d'analyses statistique R de Repenrich2 (repEnrich.R) et de RepBase (RepBase.R)

##Le dossier Visualisation

Contient les scripts pour generer les fichiers de visualisation sous IGV ou avec Circos

## Le dossier Comptages