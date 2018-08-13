################################################################
#Modification des fichiers bam pour visualisation sous IGV

#Les fichiers bam directs apres mapping ainsi que les fichiers
#filtres sont modifies

#Modification 1 : IP-Input avec bamCompare
#Modification 2 : IP/tailleLib avec bamCoverage

#Remarque : la taille des bins a ete laissee par defaut (55Kb)
################################################################


samtools index A878C21.sort.bam
samtools index A878C21.rmdup_paired_uniques.bam


#Visualisation du signal IP-Input

sample_list=(A878C17 A878C18 A878C19 A878C20)
cd ${i}
#Pour les fichiers bruts en sortie de mapping
samtools index ${i}.sort.bam
bamCompare -b1 ${i}.sort.bam -b2 /rddm1/mbertrand/Bureau/A878C21/A878C21.sort.bam --operation subtract -o ${i}BCompare.bw -p 10
#Pour les fichiers filtres (=reads uniques)
samtools index ${i}.rmdup_paired_uniques.bam
bamCompare -b1 ${i}.rmdup_paired_uniques.bam -b2 /rddm1/mbertrand/Bureau/A878C21/A878C21.rmdup_paired_uniques.bam --operation subtract -o ${i}BCompareFiltre.bw -p 10
cd ..
done

#Remarque : --operation lo2ration pour normaliser par le log2 du ration et --operation ratio pour normaliser par la taille des librairies => ces tests ont ete realises
#pour regarder si cela avait un impact sur la visualisation de nos donnees et les fichiers obtenus sont dans /analysis/IGV/regions_uniques/testBamCompare


#Visualisation du signal RPKM (IP/Taille librairie)
sample_list=(A878C17 A878C18 A878C19 A878C20)
cd ${i}
#Fichiers bruts
samtools index ${i}.sort.bam
bamCoverage -b ${i}.sort.bam  --normalizeUsing RPKM -o ${i}bCov.bw -p 10
#Fichiers filtres (=reads uniques)
samtools index ${i}.rmdup_paired_uniques.bam
bamCoverage -b ${i}.rmdup_paired_uniques.bam  --normalizeUsing RPKM -o ${i}bCovFiltre.bw -p 10
cd ..
done
