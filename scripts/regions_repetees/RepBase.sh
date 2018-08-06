#!/bin/bash


############################################################################################
#2 mars 2018
#Bertrand Mathilde

############################################################################################
#Pipeline Mise en place du protocole pour letude des regions repetees
#Alignement des reads contre RepBase en single end
#Et analyse R des donnees
############################################################################################

#1. On coupe les mates pour obtenir des tailles de 35 nt

#2. ALignement contre RepBase en single End En Very sensitive

#3. On compte combien de reads salignent sur chaque elements
#3. Analyses  : 
	#Creation dune table pour MAplot


awk -F "\t" '{if ($2 <0) absfc = -$2;else absfc = $2; 
    if (($2<=-5) && ($3>=2)) print $0 "\tblue"; 
    else if (($2>=3) && ($3>=2)) print $0 "\tred"; 
    else if ($3>=11) print $0 "\tblue"; 
    else if (($3>=3) && ($2 <= -4.5))  print $0 "\tblue"; 
    else if (($3>=3) && ($2 >= 4))  print $0 "\tred"; 
    else if (($3>=6) && ($2 <= -2))  print $0 "\tblue"; 
    else if (($3>=5) && ($2 >= 1))  print $0 "\tred"; 
    else if (($3>=8) && ($2 <= -1))  print $0 "\tblue"; 
    else if (($3>=8) && ($2 >= 1))  print $0 "\tred"; 
    else if (($3>=10) && ($2 <= -0.5))  print $0 "\tblue"; 
    else if (($3>=7) && ($2 >= 0.5))  print $0 "\tred"; 
    else print $0 "\tblack"}' MaplotRepBase.txt > MaplotRepBase.txt2
  
  R
  library(FactoMineR)


#####################################
# RepBase MAplot
#####################################
data = read.table("MaplotRepEnrich2.txt2", header = T, row.names=1)

colors = as.character(data[,3])
colors = replace(colors, colors=="red", "#d7191c")
colors = replace(colors, colors=="blue", "#2b83ba")

png("MAplot_RepEnrich2.png", width=800, height = 600)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(data[,2], data[,1], xlab="Moyenne (log2)", ylab="Fold Change (log2)", col=colors, pch = 20, cex.lab=1.5)
abline(h=0, lty=2)
abline(v=2, lty=2)
text(x=16, y=-3, label="GSAT_MM", col="#2b83ba", cex=1.5)
text(x=14, y=0.3, label="IAPEz-int", col="#2b83ba", cex=1.5)
text(x=13, y=-1.7, label="SYNREP_MM", col="#2b83ba", cex=1.5)
text(x=7, y=-4.2, label="B1_Mus1", col="#2b83ba", cex=1.5)
text(x=5.5, y=-8, label="L1_Rod", col="#2b83ba", cex=1.5)

text(x=8, y=1.5, label="B2_Mm2", col="#d7191c", cex=1.5)
text(x=6.5, y=3, label="RLTR22_Mur", col="#d7191c", cex=1.5)
#text(x=8.5, y=4, label="L1_Mur2_orf2", col="#d7191c", cex=1.5)
dev.off()


#####################################
# RepBase ACP
#####################################

data = read.table("repbase_bt2veryS_all.count", header=T, row.names=1)

res.pca = PCA(t(data[,1:4]), scale.unit=TRUE, ncp=5, graph=F)
conditions_legend = c("TALE-HD", "TALE-HDm")
colors_legend=c("#d7191c", "#2b83ba")
colors_pca= c("#d7191c", "#d7191c", "#2b83ba", "#2b83ba")


pdf("RepBase_PCA_scaleUnit.pdf", onefile = TRUE)

for (i in 1:2){
        shift = i+1
        for (j in shift:3){

                plot.PCA(res.pca, choix="ind", habillage="ind", axes=c(i,j), cex=1, col.hab=colors_pca, 
			title=paste("PCA dim",i, "/dim", j, sep=""), new.plot=$
                legend("bottomleft", legend=conditions_legend, col=colors_legend, pch=15, cex =1, border="white")

        }
}
dev.off()


    