#!/bin/bash


############################################################################################
#2 mars 2018
#Bertrand Mathilde

############################################################################################
#Pipeline Mise en place du protocole pour letude de la redistribution de H3K9me3 sur regions repetees
#Alignement des reads contre RepBase en single end
#Et analyse R des donnees
############################################################################################


#Indexation du genome pour les elements repetes : 
cd /home/raw_data/genome_ref/
mkdir genomeRepeat
cd genomeRepeat
bowtie2-build rodrep.ref rodrep
cd /home/



sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)
mate=(R1 R2)

for i in ${sample_list[*]};do
for j in ${mate[*]};do

#1. On coupe les mates pour obtenir des tailles de 35 nt
fastx_trimmer -Q33 -l 36 -i /raw_data/labo_curie/Trimme/${i}.${j}.Trimmed.fastq -o /raw_data/${i}.${j}.Trimmed35.fastq

#2. ALignement contre RepBase en single End En Very sensitive (on le fait avec les deux mates pour verifier que on a les memes resultats)
bowtie2 --very-sensitive -x /raw_data/genome_ref/genomeRepeat/rodrep -k 2 -p 4 -1 /raw_data/${i}.${j}.Trimmed35.fastq -S ${i}Rebpase_${j}.sam

#3. On compte le nbr de reads qui salignent sur chaque elements repetee : (rodrep.bed est le fichier annote repeat_maske recupere sur UCSC)
bedtools coverage -a rodrep.bed -b ${i}Rebpase_${j}.bam > ${i}_${j}RepBase.count

done
done


#Creation dune table de comptage contenant toutes les lib : 
python compile_counts.py -o repbase_R1.count A878C17_R1RepBase.count A878C18_R1RepBase.count A878C19_R1RepBase.count A878C20_R1RepBase.count A878C21_R1RepBase.count
python compile_counts.py -o repbase_R2.count A878C17_R2RepBase.count A878C18_R2RepBase.count A878C19_R2RepBase.count A878C20_R2RepBase.count A878C21_R2RepBase.count

#Analyse R des resultats : 
R
../regions_repetees/RepBase.R