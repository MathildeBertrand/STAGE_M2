#########################################################################################@
#Juillet 2018
#Tracer des violin plot pour determiner les tailles des pics EPIC

#Input : les pics de la librairie HD2 et les pics des librairies HDm mis en communs
#Pour le merge de pics HDm : concantenation, tri des chrs/start/stop et bedtools merge

#Output : un violin plot pour HD2 et un pour HDm 
###########################################################################################

library("ggplot2")

HD2=read.delim('/Users/mathildebertrand/Downloads/peak18FDR.txt', header=F,sep=" ") #Fichier de pics de HD2
HDm=read.delim('/Users/mathildebertrand/Downloads/merge_peaksHDm_e20.bed', header=F,sep="\t") #Fichier de pics de HDm

traie=(HD2$V4-HD2$V3) #Sum=387699728
Controle=c(HDm$V3-HDm$V2)#Sum=475099800

data_final_compile = data.frame(rbind(cbind("HD", as.numeric(traie)), cbind("HDm", as.numeric(Controle))))
colnames(data_final_compile) = c("DataType", "Length")
data_final_compile$Length = as.numeric(as.character(data_final_compile$Length))
data_final_compile$DataType = as.factor(data_final_compile$DataType)

data_summary <- function(x) {
  mu <- median(x)
  sigma1 <- mu-sd(x)
  sigma2 <- mu+sd(x)
  return(c(y=mu,ymin=sigma1,ymax=sigma2))
}

png("violinPlot_lengthEPIC2.png", width=600, height=600)
theme_set(
  theme_bw()
)
ggplot(data=data_final_compile, aes(x=DataType, y=Length, fill=DataType)) +  
  geom_violin(trim = FALSE) + stat_summary(fun.data=data_summary)+ 
  scale_fill_manual(values = c("#d7191c", "#2b83ba"))  + 
  labs(x = "Fréquence", y="Taille des pics EPIC") +
  theme(axis.title = element_text(size=28))  +
  theme(axis.text.y = element_text(size=20,angle=45)) +
  theme(axis.text.x = element_text(size=20))+
  theme(legend.position = "none")
dev.off()

###########################################
#Combien de pics commmuns aux 2 fichiers ?
#En commun : bedtools intersect entre les deux fichiers transformes au format bed
#Pics uniques à HD : 94039902 => 94 Mb
#Pics communs aux deux : 293659826 => 293 Mb
#Pics uniques à HDm : 181439974 =>181 Mb
###########################################
merged=read.delim('/Users/mathildebertrand/Downloads/merge.bed', header=F,sep="\t") #Fichier de pics de HDm
m=merged$V3-merged$V2

