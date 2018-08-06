################################################
#Analyse r des fichiers de sortie de RepEnrich2
################################################

library(FactoMineR)
library(RColorBrewer)
library(gplots)
library(DESeq2)

###############################################################
#Formatage des fichiers : construction dune table de comptage
#pour chacune des conditions
###############################################################

#Fichiers de donnees
HD_r1<- read.delim('/home/mbertrand/Bureau/A878C17_fraction_counts.txt', header=FALSE)
HD_r2 <- read.delim('/home/mbertrand/Bureau/A878C18_fraction_counts.txt', header=FALSE)
controle_r1 <- read.delim('/home/mbertrand/Bureau/A878C19_fraction_counts.txt', header=FALSE)
controle_r2 <- read.delim('/home/mbertrand/Bureau/A878C20_fraction_counts.txt', header=FALSE)

#Construction d'une table de comptage 
data <- data.frame(
  row.names = HD_r1[,1],
  HD_rep1 = HD_r1[,4], HD_rep2 = HD_r2[,4],
  controle_rep1 = controle_r1[,4], controle_rep2 = controle_r2[,4]
)
colData = data.frame(conditions = factor(colnames(data)))



################################################
# Strategie 1 : Analyse differentielle basique
################################################

#ACP avant DESeq
dds = DESeqDataSetFromMatrix(data,colData, formula(~conditions))
dds <- estimateSizeFactors(dds) #Estimation taille librairie
dds <- estimateDispersions(dds) #Estimation parametre de dispersion

vsd<-varianceStabilizingTransformation(dds)
data_vsd = assay(vsd)
colnames(data_vsd) = colnames(data)
hmcol = colorRampPalette(brewer.pal(5, "GnBu"))(100)

conditions_legend = c("treat", "ctrl")
colors_legend=c("#542788", "#8073ac")
colors_pca= c("#542788","#542788","#8073ac","#8073ac")

transpose = t(data_vsd)
res.pca = PCA(transpose, scale.unit=FALSE, ncp=7, graph=F)

pdf("myPCA.pdf", onefile = TRUE)
for (i in 1:6){
  shift = i+1
  for (j in shift:7){
    
    plot.PCA(res.pca, choix="ind", habillage="ind", axes=c(i,j), cex=1, col.hab=colors_pca, title=paste("PCA dim",i, "/dim", j, sep=""), new.plot=FALSE)
    legend("topright", legend=conditions_legend, col=colors_legend, pch=15, cex =1, border="white")
    
  }
}
dev.off()


######################
#DesEq
######################
data2compare = data[, c(1,2,3,4) ]

colData = data.frame(conditions = factor(c("ctrl", "ctrl", "treat", "treat")))

dds = DESeqDataSetFromMatrix(data2compare, colData, formula(~conditions))
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds = DESeq(dds)
res=results(dds)

# MAplot
#png(file="/home/mbertrand/Bureau/MAplot_DESeq2_data2compare_repEnrich2.png", height = 400, width = 800)
plotMA(dds,main="MAplot with DESeq2",ylim=c(-2,2))
#dev.off()

# Triez le fichier selon les p-values ajustée
resOrdered<-res[order(res$padj),]

# Sauvez le fichier statistique
write.table(as.data.frame(resOrdered), file="home/mbertrand/Bureau/diffana_DESeq2_data2compare_repEnrich2.txt", sep="\t")


########################################################################################
#Nouvelle strategie...
#Comme nos donnees HD1 et HD2 ne se comportent pas comme des réplicats, on change de strategie
########################################################################################


#On refais l'ACP mais en centree reduite
transpose = t(data)
res.pca = PCA(transpose, scale.unit=T, ncp=3, graph=F)

conditions_legend = c("ctrl", "treat", "public")
colors_legend=c("#542788", "#8073ac", "#b2abd2")
colors_pca= c("#542788","#542788","#8073ac","#8073ac", "#b2abd2")

pdf("myPCA.pdf", onefile = TRUE)
for (i in 1:3){
  shift = i+1
  for (j in shift:4){
    
    plot.PCA(res.pca, choix="ind", habillage="ind", axes=c(i,j), cex=1,col.hab=colors_pca, title=paste("PCA dim",i, "/dim", j, sep=""), new.plot=FALSE)
    legend("topright", legend=conditions_legend, col=colors_legend, pch=15, cex =1, border="white")
    
  }
}
dev.off()

########################################################################################
#L'ACP montre que HD1 et HD2 ne sont pas des replicats. Lanalyse precedente avec DeSeq nest donc pas envisageable
#HD1 est retiree du jeu de donnee
#On va donc comparer HD2 a la moyenne de HDm1 et HDm2
#On va juste faire un MAplot : il faut les RPKM de chaque librairie
########################################################################################

#RPKM et creation dun dataframe les contenant
HD=((HD_r2[,4]+1)/54345988)*10000000
C1=((controle_r1[,4]+1)/32999592)*10000000
C2=((controle_r2[,4]+1)/37758020)*10000000
controle=(C1+C2)/2
LogFC=log2((HD)/(controle))
moy=0.5*(HD+controle)
final=data.frame(LogFC,log2((moy)),data$type,row.names = row.names(data))

#MAplot
par(mfrow=c(1,1))
plot(((final$log2..moy..)),final$LogFC,type="p",ylab="Log2FC",xlab="Log2 de la moyenne des comptages",main="RepEnrich")

#Coloration des points
GSAT=final[grep("GSAT",row.names(final)),]
SATMIN=final[grep("SYNREP_MM",row.names(final)),]

L1=final[grep("L1",final$data.type),]#Recuperation de tous les L1
L1neg=L1
L1=L1[L1$LogFC>(0.3),]
L1=L1[L1$log2..moy..> (1),]
L1neg=L1neg[L1neg$LogFC<(-0.3),]
L1neg=L1neg[L1neg$log2..moy..> (1),]

SINE=final[grep("SINE",final$data.type),]#De tous les SINE
SINE=SINE[SINE$LogFC>(0.3),]
SINE=SINE[SINE$log2..moy..> (1),]

simple=final[grep("Simple_repeat",final$data.type),]
simpleneg=simple
simple=simple[simple$LogFC>(0.3),]
simple=simple[simple$log2..moy..> (1),]
simpleneg=simpleneg[simpleneg$LogFC< (-0.3),]
simpleneg=simpleneg[simpleneg$log2..moy..> (1),]


HAT=final[grep("hAT",final$data.type),]
HATneg=HAT
HAT=HAT[HAT$LogFC>(0.3),]
HAT=HAT[HAT$log2..moy..> (1),]
HATneg=HATneg[HATneg$LogFC<(-0.3),]
HATneg=HATneg[HATneg$log2..moy..> (1),]

tRNA=final[grep("tRNA",final$data.type),]
tRNAneg=tRNA
tRNA=tRNA[tRNA$LogFC>(0.3),]
tRNA=tRNA[tRNA$log2..moy..> (1),]
tRNAneg=tRNAneg[tRNAneg$LogFC<(-0.3),]
tRNAneg=tRNAneg[tRNAneg$log2..moy..> (1),]

ERV=final[grep("ERV",final$data.type),]
ERVneg=ERV
ERV=ERV[ERV$LogFC>(0.3),]
ERV=ERV[ERV$log2..moy..> (1),]
ERVneg=ERVneg[ERVneg$LogFC<(-0.3),]
ERVneg=ERVneg[ERVneg$log2..moy..> (1),]

points(L1$log2..moy..,L1$LogFC,col="#66a61e",pch=16) #vert a changer
points(L1neg$log2..moy..,L1neg$LogFC,col="#66a61e",pch=16)
points(simple$log2..moy..,simple$LogFC,col="#fc8d62",pch=16)#orange
points(simpleneg$log2..moy..,simpleneg$LogFC,col="#fc8d62",pch=16) #orange
points(HAT$log2..moy..,HAT$LogFC,col="#386cb0",pch=16)
points(HATneg$log2..moy..,HATneg$LogFC,col="#386cb0",pch=16) #bleu
points(tRNA$log2..moy..,tRNA$LogFC,col="#f0027f",pch=16)#rouge
points(tRNAneg$log2..moy..,tRNAneg$LogFC,col="#f0027f",pch=16)
points(GSAT$log2..moy..,GSAT$LogFC,col="#984ea3",pch=16)#vert
points(SATMIN$log2..moy..,SATMIN$LogFC,col="#984ea3",pch=16)#vert
points(ERV$log2..moy..,ERV$LogFC,col="#ffff99",pch=16)#jaune
points(ERVneg$log2..moy..,ERVneg$LogFC,col="#ffff99",pch=16)#jaune
abline(h=0,lty=2)

legend("topright", legend=c("L1", "simple_repeat","hAT","tRNA","ERV","Satellites"),col=c("#66a61e", "#fc8d62","#386cb0","#f0027f","#ffff99","#984ea3"),lty=1,lwd=3, cex=1,box.lty=0)


