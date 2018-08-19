##########################################################################
#Mathilde Bertrand
#Juillet 2018
#Y a-t-il une modification du nombre 
#de reads sur le GSAT du chr9 dans la condition traitee par rapport a la ocndition controle ?
#Si oui, peut-on la quantifier ?
##########################################################################


sample_list=(A878C17 A878C18 A878C19 A878C20 A878C21)

for i in ${sample_list[*]};do

#1 : bamCoverage (RPKM) sur toutes les libriaires avec un bs de 200
samtools index ${i}.sort.bam
bamCoverage -bs 200 --normalizeUsing RPKM -b ${i}.sort.bam -o ${i}bC_200.bw

#2. Transfromation des fichiers bw en bedgraph
./bigWigToBedGraph ${i}bC_200.bw ${i}bC_200.bedGraph

#3. On ne recupere que la region du chrs9 jusquà 34 Kb = c'est la region du GSAT 
grep (^9) ${i}bC_200.bedGraph > chr9_${i}.txt
awk '$2 <= 3040000{print $1,$2,$3,$4}' chr9_${i}.txt > GSATchr9_${i}.txt
rm chr9_${i}.txt

done

#4. Graphique avec R

R

HD_r1<- read.delim('GSATchr9_A878C17.txt', header=FALSE)
HD_r2 <- read.delim('GSATchr9_A878C18.txt', header=FALSE)
controle_r1 <- read.delim('GSATchr9_A878C19.txt', header=FALSE)
controle_r2 <- read.delim('GSATchr9_A878C20.txt', header=FALSE)

#Construction dune table de comptage 
data <- data.frame(
  start = HD_r1[,2],
  stop = HD_r1[,3],
  HD_rep1 = HD_r1[,4], 
  HD_rep2 = HD_r2[,4],
  controle_rep1 = controle_r1[,4], 
  controle_rep2 = controle_r2[,4]
)

HD2_HDm=data$HD_rep2/(0.5(data$controle_rep1+controle_rep2))
HD1_HDm=data$HD_rep1/(0.5(data$controle_rep1+controle_rep2))


plot(HD2_HDm,data$start)
lines(HD1_HDm,col='red')