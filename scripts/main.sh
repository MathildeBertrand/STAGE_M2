#!/bin/bash

###################################################
#protocole danalyse pour les donnees de ChIP-Seq
#analyse qualite des donnees
#Choix de la strategie :
#0 pour analyse globale sur genome de reference
#1 pour analyse sur les regions uniques
#2 pour analyse sur les elements repetes

#Compilation : bash main.sh
###################################################
 
#sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)
#for i in ${sample_list[*]};do
#mkdir ${i}
#./fastqc /rddm1/mbertrand/SRR6228882_1.fastq -o /rddm1/mbertrand/FASTQC/input
#fastx_trimmer -Q33 -f 10 -i /home/mbertrand/Bureau/data/nos_datas/ ${dossier_name[$j]}/${i}.fastq -o ./${i}.Trimmed.fastq
#./fastqc /home/mbertrand/Bureau/data/nos_datas/fastq_trimmer/*fastq -o /home/mbertrand/Bureau/analyse/FASTQC/nos-datas/fastqc_afterTrimmer/
#done

if [ $1 == "0" ] #On lance lanalyse globale sur le genome de reference
then 
sample_list=(A878C17.R1 A878C17.R2 A878C18.R1 A878C18.R2 A878C19.R1 A878C19.R2 A878C20.R1 A878C20.R2 A878C21.R1 A878C21.R2)
dossier_name=(A878C17 A878C18 A878C19 A878C20 A878C21)
j=0

#Trimmer le debut des sequences pour tous les fichiers
cd ../raw_data/
for i in ${sample_list[*]};do
fastx_trimmer -Q33 -f 5 -i ${i}.fastq -o ${i}.Trimmed.fastq
j=$(($j + 1))
done
#Relancer fastqc pour regarder limpact de trimmer
./Tools/fastqc /raw_data/*Trimmed.fastq -o /Analyses/FASTQC/nos-datas/fastqc_afterTrimmer/


elif [ $1 == "1" ] #On lance lanalyse sur les regions uniques
then
bash regions_uniques/filtre_strategie1.sh

elif [ $1 == "2" ] #On lance lanalyse sur les regions repetees avec RepEnrich2
then
bash regions_repetees/repEnrich.sh


else
echo "erreur : Tapez 0 pour l'analyse globale sur le genome de reference, 1 pour analyse des regions uniques et 2 pour analyses elements repetes"
fi 
