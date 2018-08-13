#!/bin/bash

##################
#27 avril 2018
#Bertrand Mathilde
#################### Protocol RepEnrich 2 #################################
#Automatisation pour l'ensemble des donnees
########################################################################

########################################################################
#Etape 1 : run the setup for RepEnrich2
########################################################################
#Le script a ete modifie de maniere a pouvoir utiliser notre version de bowtie2 sans l'option --threads
#En recuperant lannoatation via repeatmasker et via ensembl, il faut faire attention a avoir le meme nom des chrs dans les deux fichiers
#Avoir la version samtools1.3 ou plus pour que ca fonctionne correctement

python RepEnrich2_setup.py repeatmasker.txt mouse.fa dossierdereception

########################################################################
#Etape 2 : Alignement des donnees sur le genome de reference
########################################################################

#Mapping avec Bowtie2
sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)
for i in ${sample_list[*]};do
bowtie2 -q -p 6 -x /home/mbertrand/genome_ref/mouse -1 /rddm1/mbertrand/${i}.R1.Trimmed.fastq -2 /rddm1/mbertrand/${i}.R2.Trimmed.fastq -S /rddm1/mbertrand/repEnrich/${i}_mapped.sam
samtools view -bS ${i}_mapped.sam > ${i}_mapped.bam
rm ${i}_mapped.sam
done

#Lancement du script de RepEnrich2 pour distinguer les reads uniques des reads multiples
sample_list=(A878C18 A878C19 A878C20 A878C21)
for i in ${sample_list[*]};do
python ../Bureau/RepEnrich2/RepEnrich2_subset.py ${i}.bam 30 ${i} --pairedend TRUE
done

########################################################################
#Etape 3 : Run RepEnrich2 on the data => 1 semaine par librairie
########################################################################

sample_list=(A878C18 A878C19 A878C20 A878C21)
for i in ${sample_list[*]};do
python RepEnrich2.py mm10_repeatmasker.txt stepFour/ ${i} /rddm1/duvernois/STAGE_MATHILDE/RepEnrich2/results ${i}_multimap_R1.fastq --fastqfile2 ${i}_multimap_R2.fastq ${i}_unique.bam --cpus 16 --pairedend TRUE
done

########################################################################
#Etape 4 : Processing the output for RepEnrich2
########################################################################

R
repEnrich.R
