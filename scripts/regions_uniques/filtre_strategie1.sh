#!/bin/bash

############################################################################################
#2 mars 2018
#Bertrand Mathilde

############################################################################################
#Pipeline Mise en place etude des regions uniques
#Lancer : bash filtre_strategie1.sh

############################################################################################

#Etape 1 : Allignement sur le genome de reference MM10 
#Etape 2 : Les differents filtres de retrait de la redondance
#Filtre1 : on conserve uniquement les reads concordants
#Filtre2 : retrait des biais de PCR
#Filtre3 : retrait des hits multiples

#Etape 3 : On lance la detection de pics avec EPIC

#Le genome de reference et les donnees au format fastq doivent etre places dans le dossier raw_data
#Le dossier Analyses et present au meme niveau que raw_data et contient les dossiers A878C17,A878C18, A878C19, A878C20, A878C21
############################################################################################

###################
#Mapping, filtres 
###################

sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)

for i in ${sample_list[*]};do

#Alignement sur genome de reference
cd raw_data
bowtie2 --sensitive -x /genome_ref/mouse -k 2 -p 4 -1 ${i}.R1.Trimmed.fastq -2 ${i}.R2.Trimmed.fastq -S ../analysis/mapping/${i}/${i}.sam

cd ../analysis/mapping

#Application des filtres
cd ${i}
samtools sort ${i}.bam > ${i}.sort.bam #Tri des fichiers
samtools view -f 0x2 -b ${i}.sort.bam > ${i}.paired.bam #Filtre1
samtools rmdup -S ${i}.paired.bam ${i}.rmdup_paired.bam #Filtre2
samtools view -H ${i}.rmdup_paired.bam > header.sam #Filtre3
samtools view -F 4  ${i}.rmdup_paired.bam | grep -v "XS:" | cat header.sam - | samtools view -b - > ${i}.rmdup_paired_uniques.bam #Filtre3
cd ../
echo Fin de lalignement sur le genome de reference, de lapplication des filtres pour ${i}
done
cd ../
###################
#Detection de pics
###################

#Detection de pics avec macs2 pour les fichiers filtres
#macs2 callpeak -t ${i}bis.rmdup_paired_uniques.bam -c ../A878C21/A878C21bis.rmdup_paired_uniques.bam -f BAMPE -g mm -n ${i} --broad
cd Pics
mkdir EPIC
cd EPIC
#Detection de pics avec EPIC
samtools sort -n analysis/mapping/A878C21/A878C21.rmdup_paired_uniques.bam > A878C21.rmdup_paired_uniques_sort.bam
bamToBed -i A878C21.rmdup_paired_uniques_sort.bam -bedpe > A878C21.rmdup_paired_uniques_sort.bedpe
awk -F "\t" '{print "chr" $1 "\t" $2 "\t" $3 "\tchr" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' A878C21.rmdup_paired_uniques_sort.bedpe > A878C21.rmdup_paired_uniques_sortChr.bedpe

sample_list=(A878C18 A878C19 A878C20 A878C17)
for i in ${sample_list[*]};do
cd Analyses/${i}
samtools sort -n ${i}.rmdup_paired_uniques.bam> ${i}.rmdup_paired_uniques_sort.bam
bamToBed -i ${i}.rmdup_paired_uniques_sort.bam -bedpe > ${i}.rmdup_paired_uniques_sort.bedpe
awk -F "\t" '{print "chr" $1 "\t" $2 "\t" $3 "\tchr" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' ${i}.rmdup_paired_uniques_sort.bedpe > ${i}.rmdup_paired_uniques_sortChr.bedpe
epic --treatment ${i}.rmdup_paired_uniques_sortChr.bedpe --control ../A878C21.rmdup_paired_uniques_sortChr.bedpe --number-cores 20 --genome mm10 --window-size 20000 --outfile ${i}bis_unique_W20k.txt 
echo Fin EPIC  pour ${i}
cd ..
done


#Puis comptages et comparaisons des genes identifies


samtools index A878C21/A878C21bisb.rmdup_paired_uniques.bam 
sample_list=(A878C17 A878C18 A878C19 A878C20 )
for i in ${sample_list[*]};do
samtools index ${i}/${i}.rmdup_paired_uniques.bam 
echo fin index
bamCompare -b1 ${i}/${i}.rmdup_paired_uniques.bam   -b2 A878C21/A878C21.rmdup_paired_uniques.bam --scaleFactorsMethod readCount --operation subtract -p 10 -o ${i}correction.bw
done





sample_list=(A878C18 A878C19 A878C20 A878C17)
for i in ${sample_list[*]};do
cd ${i}
bamCompare -b1 ${i}bis.sort.bam -b2 /rddm1/mbertrand/Bureau/A878C21/A878C21bis.sort.bam --operation subtract -bs 5 -o ${i}BComparebs5.bw -p 10
bamCoverage -b ${i}bis.sort.bam  -bs 5 --normalizeUsing RPKM -o ${i}bCovbs5.bw -p 10
bamCoverage -b ${i}bis.sort.bam  -bs 10 --normalizeUsing RPKM -o ${i}bCovbs10.bw -p 10

bamCompare -b1 ${i}bisb.rmdup_paired_uniques.bam -b2 /rddm1/mbertrand/Bureau/A878C21/A878C21bisb.rmdup_paired_uniques.bam --operation subtract -bs 5 -o ${i}BComparebs5Filtre.bw -p 10
bamCoverage -b ${i}bisb.rmdup_paired_uniques.bam  -bs 5 --normalizeUsing RPKM -o ${i}bCovbs5Filtre.bw -p 10
bamCoverage -b ${i}bisb.rmdup_paired_uniques.bam  -bs 10 --normalizeUsing RPKM -o ${i}bCovbs10Filtre.bw -p 10

cd ..
done
