##############################################################################################################################
#Mathilde Bertrand
#Juin/Juillet 2018

#Analyse R des resultats de RepBase : ACP et MAplot
#Ce que l'on veut faire : on veut determiner si on a un enrichissement ou une perte de la methylation sur les elements repetes
#On veut comparer ces resultats avec ceux que lon obtient avec RepEnrich2

#Input : le fichier de comptages du nombre de reads sur chaque element repetes obtenus en sortie de RepBase avec les lib 17, 18, 19 et 20
#Remarque : il est plus difficile de colorier les points ici qu'avec RepEnrich2 car les elements repetes (meme si sont les memes)
#nont pas le meme noms et il est plus dur dautomatiser la coloration
##############################################################################################################################

library(FactoMineR)


count=read.delim('/Users/mathildebertrand/Downloads/repbase_R1.count', header=T)

######################################################
#ACP pour regarder le comportement des librairies :
# en centree reduite
######################################################

transpose = t(count)
res.pca = PCA(transpose, scale.unit=T, ncp=3, graph=F)

conditions_legend = c("ctrl", "treat")
colors_legend=c("blue", "red")
colors_pca= c("red","red","blue","blue")

pdf("myPCA.pdf", onefile = TRUE)
for (i in 1:3){
  shift = i+1
  for (j in shift:4){
    plot.PCA(res.pca, choix="ind", habillage="ind", axes=c(i,j), cex=1,col.hab=colors_pca, title=paste("PCA dim",i, "/dim", j, sep=""), new.plot=FALSE)
    legend("topright", legend=conditions_legend, col=colors_legend, pch=15, cex =1, border="white")
  }
}
dev.off()


###########################################################################################
#HD1 semble se comporter comme un controle. On lenleve du jeu de donnees
#On ne peut plus faire danalyse statistique, on se limite donc a la realisation dun MAplot
###########################################################################################

HD=((count$A878C18+1)/54345988)*10000000
C1=((count$A878C19+1)/32999592)*10000000
C2=((count$A878C20+1)/37758020)*10000000
controle=(C1+C2)/2
LogFC=log2((HD)/(controle))
moy=0.5*(HD+controle)
final=data.frame(LogFC,log2((moy)),row.names = count$element)
plot(((final$log2..moy..)),final$LogFC,type="p",ylab="Log2FC",xlab="Log2 de la moyenne des comptages",main="RepBase")
abline(h=0,lty=2)

#---------------------
#Coloriage
#---------------------

#Recuperation de tous les L1
L1=final[grep("^L1",row.names(final)),]
L1neg=L1
L1=L1[L1$LogFC>(0.3),]
L1=L1[L1$log2..moy..> (1),]
L1neg=L1neg[L1neg$LogFC<(-0.3),]
L1neg=L1neg[L1neg$log2..moy..> (1),]

#Recuperation des satellites
GSAT=final[grep("GSAT",row.names(final)),]
SATMIN=final[grep("SATMIN",row.names(final)),]

#tRNA
tRNA=final[grep("tRNA",row.names(final)),]
tRNAneg=tRNA
tRNA=tRNA[tRNA$LogFC>(0.3),]
tRNA=tRNA[tRNA$log2..moy..> (1),]
tRNAneg=tRNAneg[tRNAneg$LogFC<(-0.3),]
tRNAneg=tRNAneg[tRNAneg$log2..moy..> (1),]

#ERV
ERV=final[grep("ERV",row.names(final)),]
ERVneg=ERV
ERV=ERV[ERV$LogFC>(0.3),]
ERV=ERV[ERV$log2..moy..> (1),]
ERVneg=ERVneg[ERVneg$LogFC<(-0.3),]
ERVneg=ERVneg[ERVneg$log2..moy..> (1),]

points(2.615493,0.6789171,col="#fc8d62",pch=16)
points(2.615493,0.6789171,col="#fc8d62",pch=16)
points(L1$log2..moy..,L1$LogFC,col="#66a61e",pch=16) #vert a changer
points(L1neg$log2..moy..,L1neg$LogFC,col="#66a61e",pch=16)
points(GSAT$log2..moy..,GSAT$LogFC,col="#984ea3",pch=16)#vert
points(SATMIN$log2..moy..,SATMIN$LogFC,col="#984ea3",pch=16)#vert
points(ERV$log2..moy..,ERV$LogFC,col="#ffff99",pch=16)#jaune
points(ERVneg$log2..moy..,ERVneg$LogFC,col="#ffff99",pch=16)#jaune

legend("topright", legend=c("L1", "simple_repeat","hAT","tRNA","ERV","Satellites"),col=c("#66a61e", "#fc8d62","#386cb0","#f0027f","#ffff99","#984ea3"),lty=1,lwd=3, cex=1,box.lty=0)






