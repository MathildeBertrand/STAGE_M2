########################################################
#Fin juin 2018
#Mathilde Bertrand
#Analyse du genome decoupe en fenetres
#Et compter le nombre de reads dans chacune des fenêtres
########################################################

#On determine la longueur et le nombre de fenetre que lon va utiliser : 
python /home/duvernois/SOFTS/UTILS/get_length_fasta.py -f Mus_musculus.GRCm38.dna.primary_assembly.fa -o Mus_musculus.GRCm38.dna.primary_assembly.count
# genome size : 2 730 871 774
# pour 50 000 fenetres ==> fenetres d'environ 55Kb

#################################################################
# utilisation des bedtools pour creer les fenetres
#Fenetres de 55, 15, 15 et 5 Kb
# => Pour regarder si la taille des fenetres change les resultats
#################################################################

sed '1d' Mus_musculus.GRCm38.dna.primary_assembly.count | awk -F "\t" '{print $1 "\t" $2}' > Mus_musculus.GRCm38.dna.primary_assembly.genome

sample_list=(55 25 15 5)
for i in ${sample_list[*]};do
bedtools makewindows -g Mus_musculus.GRCm38.dna.primary_assembly.genome -w ${i} > Mus_musculus.GRCm38.dna.primary_assembly.win${i}K.bed
done

awk -F "\t" '{print $1 "\tmm10\twin55K\t" $2 "\t" $3 "\t.\t+\t.\tID=win55Kb_" NR}' Mus_musculus.GRCm38.dna.primary_assembly.win55K.bed | awk -F "\t" '{if ($4 ==0) print $1 "\t" $2 "\t" $3 "\t1\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9; else print $0}' > Mus_musculus.GRCm38.dna.primary_assembly.win55K.gff
awk -F "\t" '{print $1 "\tmm10\twin25K\t" $2 "\t" $3 "\t.\t+\t.\tID=win25Kb_" NR}' Mus_musculus.GRCm38.dna.primary_assembly.win25K.bed | awk -F "\t" '{if ($4 ==0) print $1 "\t" $2 "\t" $3 "\t1\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9; else print $0}' > Mus_musculus.GRCm38.dna.primary_assembly.win25K.gff
awk -F "\t" '{print $1 "\tmm10\twin15K\t" $2 "\t" $3 "\t.\t+\t.\tID=win15Kb_" NR}' Mus_musculus.GRCm38.dna.primary_assembly.win15K.bed | awk -F "\t" '{if ($4 ==0) print $1 "\t" $2 "\t" $3 "\t1\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9; else print $0}' > Mus_musculus.GRCm38.dna.primary_assembly.win15K.gff
awk -F "\t" '{print $1 "\tmm10\twin5K\t" $2 "\t" $3 "\t.\t+\t.\tID=win5Kb_" NR}' Mus_musculus.GRCm38.dna.primary_assembly.win5K.bed | awk -F "\t" '{if ($4 ==0) print $1 "\t" $2 "\t" $3 "\t1\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9; else print $0}' > Mus_musculus.GRCm38.dna.primary_assembly.win5K.gff


###############################################
#Comptage
###############################################
# creation d'un bam "single-end"  : 
#sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)
#for i in ${sample_list[*]};do
#bamToBed -i A878C${i}.rmdup_paired_uniques.bam | awk -F "\t" '{print $1 "\t" $2 "\t" $3 "\tread_" NR "\t" $5 "\t" $6}' | bedToBam -i stdin -g /home/duvernois/PROJECTS/STAGE_MATHILDE/cut_signal_windows/Mus_musculus.GRCm38.dna.primary_assembly.genome > A878C${i}bis.rmdup_paired_uniques.SE.bam
#done
#=>
# comptage avec htseq count => PROBLEME

################################################################################################################
# Du coup, comptages avec bedtools coverage (version 2.21) qui pour chaque fenetre calcule la couverture de reads 
################################################################################################################

sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)
for i in ${sample_list[*]};do
bedtools coverage -b Mus_musculus.GRCm38.dna.primary_assembly.win55K.bed -abam ${i}.rmdup_paired_uniques.bam > ${i}.rmdup_paired_uniques_win55K.count
bedtools coverage -b Mus_musculus.GRCm38.dna.primary_assembly.win25K.bed -abam ${i}.rmdup_paired_uniques.bam > ${i}.rmdup_paired_uniques_win25K.count
bedtools coverage -b Mus_musculus.GRCm38.dna.primary_assembly.win15K.bed -abam ${i}.rmdup_paired_uniques.bam > ${i}.rmdup_paired_uniques_win15K.count
bedtools coverage -b Mus_musculus.GRCm38.dna.primary_assembly.win5K.bed -abam ${i}.rmdup_paired_uniques.bam > ${i}.rmdup_paired_uniques_win5K.count
done

############################################################
#Creation dun fichier de comptage pour toutes les librairies
#A chaque tour de boucle, on change la taille des fenetres
############################################################
sample_list=(55 25 15 5)
for i in ${sample_list[*]};do
python compile_counts.py -o win${i}K_IP.count A878C17.rmdup_paired_uniques_win${i}K.count A878C18.rmdup_paired_uniques_win${i}K.count A878C19.rmdup_paired_uniques_win${i}K.count A878C20.rmdup_paired_uniques_win${i}K.count
done

############################################
#Analyse des resultats avec R : RPKM et ACP
############################################
R
library(FactoMineR)

data = read.table("55/win55K_IP.count", header=T, row.names=1)

data <- data.frame(
  row.names=row.names(data),
  HD1=data$A878C17bis.rmdup,
  HD2 = data$A878C18bis.rmdup,
  HDm1=data$A878C19bis.rmdup,
  HDm2=data$A878C20bis.rmdup
)


#ACP
res.pca = PCA(t(data), scale.unit=TRUE, ncp=5, graph=F)
conditions_legend = c("TALE-HD", "TALE-HDm")
colors_legend=c("#d7191c" ,"#2b83ba")
colors_pca= c("#d7191c", "#d7191c", "#2b83ba", "#2b83ba")
pdf("55/win55K_IP_PCA_scaleUnit_all.pdf", onefile = TRUE) 

for (i in 1:2){ 
        shift = i+1
        for (j in shift:3){
                plot.PCA(res.pca, choix="ind", habillage="ind", axes=c(i,j), cex=1, col.hab=colors_pca, title=paste("PCA dim",i, "/dim", j, sep=""), new.plot=FALSE)
                legend("bottomleft", legend=conditions_legend, col=colors_legend, pch=15, cex =1, border="white")
        }
}
dev.off()


#Creation table de comptage RPKM pour MAplot => permet de verifier la coherence des resultats
HD2=(data$A878C18bisb.rmdup/54345988)*10000000
HDm1=(data$A878C19bisb.rmdup/32999592)*10000000
HDm2=(data$A878C20bisb.rmdup/37758020)*10000000
HDm=0.5*(HDm1+HDm2)
LogFC=log2(HD2/HDm)
moy=0.5µ(HD2+HDm)
moy=0.5*(HD2+HDm)
final=data.frame(LogFC,log2(moy),row.names=row.names(data))
write.table(as.data.frame(final), file="Compatges15.txt", sep="\t") 










########################################################################################################
#Test dune nouvelle methode pour decouper le genome en fenetres
#Utilisation de bamCoverage pour calculer des fentres en rpkm => na pas donne des resultats concluants
#Car ne decoupe pas les fenetres de la meme maniere dun echantillon a un autre
########################################################################################################


#Test de fenetres de 5,55,15 et 25 kb :
sample_list=(A878C17 A878C18 A878C19 A878C20 )
for i in ${sample_list[*]};do
cd ${i}
bamCoverage -b ${i}bisb.rmdup_paired_uniques.bam  -bs 5 --normalizeUsing RPKM -o ${i}bCovbs5Filtre.bw -p 10
bamCoverage -b ${i}bisb.rmdup_paired_uniques.bam  -bs 55 --normalizeUsing RPKM -o Fenetres/${i}bCovbs55Filtre.bw -p 10
bamCoverage -b ${i}bisb.rmdup_paired_uniques.bam  -bs 15 --normalizeUsing RPKM -o Fenetres/${i}bCovbs15Filtre.bw -p 10
bamCoverage -b ${i}bisb.rmdup_paired_uniques.bam  -bs 25 --normalizeUsing RPKM -o Fenetres/${i}bCovbs25Filtre.bw -p 10
cd ..
done


#Conversion des fichiers en bedgraph
sample_list=(A878C17 A878C18 A878C19 A878C20)

./bigWigToBedGraph Fenetres/${i}bCovbs5Filtre.bw Fenetres/${i}bs5.bedgraph
./bigWigToBedGraph Fenetres/${i}bCovbs55Filtre.bw Fenetres/${i}bs55.bedgraph
./bigWigToBedGraph Fenetres/${i}bCovbs15Filtre.bw Fenetres/${i}bs15.bedgraph
./bigWigToBedGraph Fenetres/${i}bCovbs25Filtre.bw Fenetres/${i}bs25.bedgraph


#Creation de la table de comptage et ajout de lentete
awk -F "\t" '{print "ID"NR "\t" $1 "\t" $2 "\t" $3 "\t" $4}' HD2.bedgraph > HD22.bedgraph

R

HD_r1<- read.delim('HD1.bedgraph', header=FALSE)
HD_r2 <- read.delim('HD22.bedgraph', header=FALSE)
controle_r1 <- read.delim('HDm1.bedgraph', header=FALSE)
controle_r2 <- read.delim('HDm2.bedgraph', header=FALSE)


#Construction de la  table de comptage
data <- data.frame(
row.names = HD_r2[,1],
  HD1 = HD_r1[,4], HD2 = HD_r2[,5],
  HDm1 = controle_r1[,4], HDm2 = controle_r2[,4],
)

library(FactoMineR)

#data = read.table("ACPbamCoverage.txt", header=T, row.names=1)

res.pca = PCA(t(data), scale.unit=TRUE, ncp=5, graph=F)
conditions_legend = c("HD","HDm")
colors_legend=c("#d7191c" ,"#2b83ba")
colors_pca= c("#d7191c", "#d7191c", "#2b83ba", "#2b83ba")

pdf("BamCoverage_PCA.pdf", onefile = TRUE) 

for (i in 1:2){ 
        shift = i+1
        for (j in shift:3){
                
                plot.PCA(res.pca, choix="ind", habillage="ind", axes=c(i,j), cex=1, col.hab=colors_pca, title=paste("PCA dim",i, "/dim", j, sep=""), new.plot=FALSE)
                legend("bottomleft", legend=conditions_legend, col=colors_legend, pch=15, cex =1, border="white")
                
        }
}
dev.off()
