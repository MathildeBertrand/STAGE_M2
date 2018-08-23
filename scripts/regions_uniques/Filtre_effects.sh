######################################################################################################
#Mars 2018
#Quantifier le nombre de reads que l'on perd apres application de differents filtres de la redondance
#Pour chaque filtre, le nombre de reads present dans le fichier est affiche dans la konsole

#Outils : samtools, version 1.3.1
#Lancer : se placer dans le repertoire ou sont stockes les fichiers et tapper : bash Filtre_effects.sh
######################################################################################################



sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)

for i in ${sample_list[*]};do

echo Debut de calcul pour la librairie ${i}  :

echo Pour le fichier ${i}.sort.bam :
samtools view -F 0x904 ${i}.sort.bam | sort | uniq | wc -l #-F 0x904 : retire les reads non alignes et les reads qui sont des alignements secondaires

echo Pour le fichier ${i}.paired.bam : #Nombre de reads restant apres le premier filtre
samtools view -F 0x904 ${i}.paired.bam | sort | uniq | wc -l

echo Pour le fichier ${i}.rmdup_paired.bam : #Nombre de reads restant apres le deuxieme filtre
samtools view -F 0x904 ${i}.rmdup_paired.bam | sort | uniq | wc -l

echo Pour le fichier ${i}.rmdup_paired_uniques.bam : #Nombre de reads restant apres le troisieme filtre
samtools view -F 0x904 ${i}.rmdup_paired_uniques.bam | sort | uniq | wc -l

done