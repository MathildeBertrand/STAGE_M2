#Tracer des Manhnattan plot a partir des pics EPIC pour chaque librairie
#Le but est de regarder la repartition des pics pour enlever les moins bons et ne garder que les meilleurs


library(stringr)
library(qqman)

data= read.table("A878C18bis_unique_W20k.txt", header=T) #Fichier de pics de la lib 18
data$BP = data[,3] - data[,2]
data$SNP=paste(data$Chromosome, data$BP, sep="_")

t=str_replace(data$Chromosome, "chr","")
t2=str_replace(t, "X","20")
data$CHR=as.numeric(t2)

png("manhattanPlot_lib18EPIC.png", height=600, width=1000)
par(mar=c(5.1, 7, 4.1, 2.1))
manhattan(data, ylim=c(0,400), suggestiveline = F, genomewideline = F, 
	col=c("powderblue", "gold1"), cex=2, 
	chrlabs=c(1:19, "X"), cex.lab=2, cex.axis=1.5)
abline(h=10, lty=2)
dev.off()




data= read.table("A878C19bis_unique_W20k.txt", header=T) #Fichier de pics de la lib 19
data$BP = data[,3] - data[,2]
data$SNP=paste(data$Chromosome, data$BP, sep="_")

t=str_replace(data$Chromosome, "chr","")
t2=str_replace(t, "X","20")
data$CHR=as.numeric(t2)

png("manhattanPlot_lib19EPIC.png", height=600, width=1000)
par(mar=c(5.1, 7, 4.1, 2.1))
manhattan(data, ylim=c(0,400), suggestiveline = F, genomewideline = F, 
	col=c("powderblue", "gold1"), cex=2, 
	chrlabs=c(1:19, "X"), cex.lab=2, cex.axis=1.5)
abline(h=10, lty=2)
dev.off()





data= read.table("A878C20bis_unique_W20k.txt", header=T) #Fichier de pics de la lib 20
data$BP = data[,3] - data[,2]
data$SNP=paste(data$Chromosome, data$BP, sep="_")

t=str_replace(data$Chromosome, "chr","")
t2=str_replace(t, "X","20")
data$CHR=as.numeric(t2)

png("manhattanPlot_lib20EPIC.png", height=600, width=1000)
par(mar=c(5.1, 7, 4.1, 2.1))
manhattan(data, ylim=c(0,400), suggestiveline = F, genomewideline = F, 
	col=c("powderblue", "gold1"), cex=2, 
	chrlabs=c(1:19, "X"), cex.lab=2, cex.axis=1.5)
abline(h=10, lty=2)
dev.off()
