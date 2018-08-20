#!/bin/bash

#########################################
#Une fois lannotation des genes realisee, 
#On recupere les genes qui nous interesse
##########################################


grep -E "EID|DIE|OA|OB" outfile.txt > geneOfInterest.txt #On recupere les lignes qui sont marquees DIE ou EID
#outfile.txt est le fichier de sortie du script annotate_domains_v3.0.py

awk '{print $6}' geneOfInterest.txt > ListeGene.txt #On extrait le nom des genes marques DIE ou EID

#La liste de genes est ensuite copiee dans Venny pour comparer les genes identifies dans les differentes conditions
