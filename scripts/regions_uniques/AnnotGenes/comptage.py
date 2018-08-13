#!/bin/bash

#Mathilde Bertrand


#Combien de reads on a sur les genes ?
#A partir de la sortie de HTSEQcount qui a compte pour chaque pic cmbien il y avait de reads dessus :
#on prend le fichier de sortie de annotate domain, on recupere tous les pics qui sont directement annotes sur des genes
#n compare ces deux fichiers de maniere a recuperer pour chacun des pics annotes sur des genes combien il y a de reads dessus


grep -E '(DIE|EID)' outfileHD.txt > O.txt #on recupere tous ce qui est annote directement sur des genes
awk '{print $1}' O.txt > PicsGENes.txt #On recupere uniquement le noms de pics

for i in `cat PicsGENes.txt | awk {'print $1'}`;do
grep ${i} HTSEQcountGene/NurHD/A878C18HD.txt >> inter.txt
done
#La somme des comptages de toutes les lignes nous donne le nombre de reads qui sont sur des genes
awk '{print (total +=$2)}' inter.txt

#On fait le menage
rm PicsGENes.txt
rm O.txt
rm inter.txt

#Pour le comptage des reads en intergenique
grep -E '(DIE|EID)' outfileC.txt > O.txt #on recupere tous ce qui est annote directement sur des genes
awk '{print $1}' O.txt > PicsGENes.txt #On recupere uniquement le noms de pics

for i in `cat PicsGENes.txt | awk {'print $1'}`;do
grep ${i} HTSEQcountGene/NurC/A878C20C.txt >> inter.txt
done
#La somme des comptages de toutes les lignes nous donne le nombre de reads qui sont sur des genes
awk '{print (total +=$2)}' inter.txt

#On fait le menage
rm PicsGENes.txt
rm O.txt
rm inter.txt
