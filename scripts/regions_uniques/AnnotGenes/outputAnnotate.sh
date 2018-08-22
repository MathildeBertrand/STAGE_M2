#!/bin/bash

#########################################
#Une fois lannotation des genes realisee, 
#On recupere les genes qui nous interesse
##########################################


grep -E "EID|DIE|OA|OB" outfile.txt > geneOfInterest.txt #On recupere les lignes qui sont marquees DIE ou EID
#outfile.txt est le fichier de sortie du script annotate_domains_v3.0.py

awk '{print $6}' geneOfInterest.txt > ListeGene.txt #On extrait le nom des genes marques DIE ou EID
awk '{print $1}' picsOfInterest.txt > Listepics.txt #On extrait le nom des pics qui sont annotes sur des genes

#La liste de genes est des pics sont ensuite copiees dans Venny pour comparer les genes identifies dans les differentes conditions (ou le nb de pics dans chaque cond)
