#################################################################################################################
#Creation des fichiers d'Input necessaire a la realisation des circos : 

#1.  les fichiers de signaux ont ete generes au format bw grace a bamCompare (operation subtract et bin de 400)
#puis ont ete convertit en ou bedgraphe 
#=> sont trop denses pour quil y ait un rendu sur le graphique
#Du coup on filtre en prenant les posititions >= 10 reads et en fixant à 100 la limite supérieure.
#et on ajoute les bonnes couleurs

#2. On ajoute les bonnes couleurs pour les fichiers de pics 
#3. On cree le circos
################################################################################################################@


#Pour la librairie 17 (HD1), on la represente en rouge
awk -F "\t" '{if ($4 >=10) print $0 "\tfill_color=red"}' DATA_init/A878C17_400bC.bedGraph | awk -F "\t" '{if ($4 >= 100) print $1 "\t" $2 "\t" $3 "\t100"; else print $0}' > A878C17_400bC_10-100.txt

#La librairie 18 (HD2) sera en rouge
awk -F "\t" '{if ($4 >=10) print $0 "\tfill_color=orange"}' DATA_init/A878C18_400bC..bedGraph | awk -F "\t" '{if ($4 >= 100) print "mm" $1 "\t" $2 "\t" $3 "\t100"; else print "mm" $0}' > A878C18_400bC_10-100.txt

#HDm1 et HDm2 en bleu => ce sont les deux controles et on a vu par les ACP que les comportements etaient similaires
awk -F "\t" '{if ($4 >=10) print $0 "\tfill_color=blue"}' DATA_init/A878C19_400bC.bedGraph | awk -F "\t" '{if ($4 >= 100) print "mm" $1 "\t" $2 "\t" $3 "\t100"; else print "mm" $0}' > positive19_test.txt
awk -F "\t" '{if ($4 >=10) print $0 "\tfill_color=blue"}' DATA_init/A878C20_400bCbedGraph | awk -F "\t" '{if ($4 >= 100) print "mm" $1 "\t" $2 "\t" $3 "\t100"; else print "mm" $0}' > positive20_test.txt


#Pour les pics de macs2 : On donne les bonnes couleurs et on ajoute mm en debut pour que se positionne correctement sur les bons chrs
#awk -F "\t" '{print "mm" $1 "\t" $2 "\t" $3 "\t1\tfill_color=red"}' correctA878C17_peaks.broadPeak > correctA878C17_peaks.txt
#awk -F "\t" '{print "mm" $1 "\t" $2 "\t" $3 "\t1\tfill_color=red"}' correctA878C18_peaks.broadPeak > correctA878C18_peaks.txt
#awk -F "\t" '{print "mm" $1 "\t" $2 "\t" $3 "\t1\tfill_color=blue"}' correctA878C19_peaks.broadPeak > correctA878C19_peaks.txt
#awk -F "\t" '{print "mm" $1 "\t" $2 "\t" $3 "\t1\tfill_color=blue"}' correctA878C20_peaks.broadPeak > correctA878C20_peaks.txt


#Pour les pics EPIC : 
A878C17bis_unique_W20k_nofilter.txt
A878C18bis_unique_W20k_nofilter.txt
A878C19bis_unique_W20k_nofilter.txt
A878C20bis_unique_W20k_nofilter.txt


#On lance le circos
~/Bureau/TOOLS/circos-0.69/bin/circos --conf circos.conf

