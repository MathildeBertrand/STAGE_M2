#!/bin/bash

############################################################################################
#2 mars 2018
#Bertrand Mathilde

############################################################################################
#Logiciels necessaires : 
#bowtie2
#samtools
#EPIC
############################################################################################

############################################################################################
#Pipeline Mise en place de lanalyse des regions uniques
#Lancer :
#Via main.sh ou aussi par : bash filtre_strategie1.sh
############################################################################################

#Etape 1 : Allignement sur le genome de reference mm10 
#Etape 2 : Les differents filtres de retrait de la redondance
#Filtre1 : on conserve uniquement les reads concordants
#Filtre2 : retrait des biais de PCR
#Filtre3 : retrait des hits multiples
#Les fichiers generes seront dans  : /analysis/mapping/$data$i/

#Etape 3 : On lance la detection de pics avec EPIC (les parametres de la detection avec MACS2 ont ete conserves dans les commentaires)
#Les fichiers de sortie EPIC sont ranges dans analysis/Pics/EPIC/output_EPIC

#Le genome de reference et les donnees au format fastq doivent etre places dans le dossier raw_data
#Le dossier Analyses et present au meme niveau que raw_data et contient les dossiers A878C17,A878C18, A878C19, A878C20, A878C21
############################################################################################

#Creation des dossiers de rangements de la detection de pics
cd ../../analysis/
mkdir Pics
mkdir mapping
cd Pics/
mkdir EPIC

###################
#Mapping, filtres 
###################

sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)

for i in ${sample_list[*]};do

#Creation des dossiers de rangement de donnees de mapping et de filtre : 
if [ $i == "A878C17" ] 
then
dossier='data17'
mkdir  /analysis/mapping/data17/
#echo $i

elif [ $i == "A878C18" ] 
then
dossier='data18'
mkdir  /analysis/mapping/data18/

elif [ $i == "A878C19" ] 
then
dossier='data19'
mkdir  /analysis/mapping/data19/

elif [ $i == "A878C20" ] 
then
dossier='data20'
mkdir  /analysis/mapping/data20/

else [ $i == "A878C21" ] 
dossier='data21'
mkdir  /analysis/mapping/data21/
fi

echo Debut de lalignement sur le genome de reference et de lapplication des filtres pour ${i}

bowtie2 --sensitive -x /raw_data/genome_ref/mouse -k 2 -p 4 -1 /raw_data/labo_curie/Trimme/${i}.R1.Trimmed.fastq -2 /raw_data/labo_curie/Trimme/${i}.R2.Trimmed.fastq -S /analysis/mapping/$dossier/${i}.sam #Alignement sur genome de reference
samtools view -Sb  /analysis/mapping/$dossier/${i}.sam  > /analysis/mapping/$dossier/${i}.bam #Conversion sam en bam 

#Application des filtres 
cd ../../mapping/$dossier
samtools sort ${i}.bam > ${i}.sort.bam #Tri des fichiers

mkdir Filtres #Dossier qui contiendra les fichiers de mapping filtres
cd Filtres
samtools view -f 0x2 -b ${i}.sort.bam > ${i}.paired.bam #Filtre1
samtools rmdup -S ${i}.paired.bam ${i}.rmdup_paired.bam #Filtre2
samtools view -H ${i}.rmdup_paired.bam > header.sam #Filtre3
samtools view -F 4  ${i}.rmdup_paired.bam | grep -v "XS:" | cat header.sam - | samtools view -b - > ${i}.rmdup_paired_uniques.bam #Filtre3

echo Fin de lalignement sur le genome de reference et de lapplication des filtres pour ${i}


#########################################################
#Pre-traitements pour la detection de pics EPIC
#########################################################

echo Pre-traitements pour la detection des pics de ${i}
samtools sort -n ${i}.rmdup_paired_uniques.bam> ${i}.rmdup_paired_uniques_sort.bam #Sort des fichiers
bamToBed -i ${i}.rmdup_paired_uniques_sort.bam -bedpe > ${i}.rmdup_paired_uniques_sort.bedpe #Conversion bam en bed
awk -F "\t" '{print "chr" $1 "\t" $2 "\t" $3 "\tchr" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' ${i}.rmdup_paired_uniques_sort.bedpe > ../../../EPIC/${i}.rmdup_paired_uniques_sortChr.bedpe
echo Fin des pre-traitements pour la detection des pics de ${i}

cd ../../../../ #On se replace dans le dossier home
done



###################
#Detection de pics
###################

#Detection de pics avec macs2 pour les fichiers filtres : 
#macs2 callpeak -t ${i}bis.rmdup_paired_uniques.bam -c ../A878C21/A878C21bis.rmdup_paired_uniques.bam -f BAMPE -g mm -n ${i} --broad

#Detection de pics avec EPIC : 

cd analysis/pics/EPIC
sample_list=(A878C17 A878C18 A878C19 A878C20)

for i in ${sample_list[*]};do
echo Debut EPIC pour ${i}
epic --treatment ${i}.rmdup_paired_uniques_sortChr.bedpe --control A878C21.rmdup_paired_uniques_sortChr.bedpe --number-cores 20 --genome mm10 --window-size 20000 --outfile ${i}bis_unique_W20k.txt 
echo Fin EPIC  pour ${i}
done

#Rangements des fichiers de la detection de pics EPIC
mkdir Input #les fichiers pour lancer EPIC
mkdir output_EPIC #Les fichiers de sortie EPIC
mv *rmdup_paired_uniques_sortChr.bedpe Input/
mv *bis_unique_W20k.txt output_EPIC/



