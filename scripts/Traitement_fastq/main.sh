#!/bin/bash

#########################################################################################################################################################
#Mathilde Bertrand 2018

#Protocole danalyse pour les donnees de ChIP-Seq qui permet de generer les fichiers necessaires aux analyses R a partir des fichiers fastq
#de maniere a obtenir les resultats presentes dans le rapport de stage
#3 strategies differentes peuvent etre lancee a partir de ce script


#Choix de la strategie :
#0  : pour analyse globale sur genome de reference a savoir lanalyse qualite des sequences 
#En tapant 0, il faut que le programme FASTQC soit range dans le dossier Tools
#les fichiers produits sont des fichiers de sortie FASTQC nommes ".Trimmed.fastq" car les 5 premieres paires de bases sont retirees
#et ranges dans le dossier /Analyses/FASTQC/nos-datas/fastqc_afterTrimmer/


#1  : pour analyse sur les regions uniques :  alignement sur genome de reference, filtre de la redondance, detection de pics 
#En tapant 1, on fait appel au script "filtre_strategie1.sh" qui est range dans le dossier regions_uniques

#2  : pour analyse sur les elements repetes : RepEnrich2 sera lance ici (RepBase.sh est disponible dans le meme dossier mais comme donne les memes
#resultats que RepEnrich2, jai fait le choix de lancer RepEnrich2 par defaut).

#Compilation : 
#bash main.sh 0
#bash main.sh 1
#bash main.sh 2

#Pour que les differentes strategies soient lancees, il faut que les sequences au format FASTQ soit rangees dans le dossier raw_data
#########################################################################################################################################################
 

#############################################
#Analyse globale sur le genome de reference
#############################################

if [ $1 == "0" ] 
then 
sample_list=(A878C17.R1 A878C17.R2 A878C18.R1 A878C18.R2 A878C19.R1 A878C19.R2 A878C20.R1 A878C20.R2 A878C21.R1 A878C21.R2)
dossier_name=(A878C17 A878C18 A878C19 A878C20 A878C21)
j=0

cd ../../Tools

#Analyse qualite avec les sequences completes : 

./fastqc /raw_data/labo_curie/fastq/${i}.fastq -o /analysis/FASTQC/fastqc_beforeTrimmer/

#Trimmer le debut des sequences pour tous les fichiers
#cd ../raw_data/
for i in ${sample_list[*]};do
fastx_trimmer -Q33 -f 5 -i /raw_data/labo_curie/fastq/${i}.fastq -o /raw_data/${i}.Trimmed.fastq
j=$(($j + 1))

#Relancer fastqc pour regarder limpact de trimmer
./fastqc /raw_data/labo_curie/fastq/${i}.Trimmed.fastq -o /analysis/FASTQC/nos-datas/fastqc_afterTrimmer/
done

#Rangement : 
cd ../raw_data/labo_curie
mkdir Trimme #Creation dun dossier de rangement des fastq trimes
mv *${i}.Trimmed.fastq Trimme/


#################################
#Analyse sur les regions uniques
#################################

elif [ $1 == "1" ] 
then
bash filtre_strategie1.sh

###################################
#Analyse sur les regions repetees
###################################

elif [ $1 == "2" ] 
then
bash repEnrich.sh
#bash RepBase.sh

###################################
########Message derreur############
###################################
else
echo "erreur : Tapez 0 pour l'analyse globale sur le genome de reference, 1 pour analyse des regions uniques et 2 pour analyses elements repetes"
fi 
