#Mathilde Bertrand
#avril 2018
#Enlever le mot "chr" devant chaque chrs pour la compatibilite des fichiers

sed -e 's/chr//g' repeatAnnotation.gtf/repeatAnnotation.gtf > repeat.gtf
