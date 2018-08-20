#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Mathilde Bertrand
7 mai 2018

 conversion dun fichier de pics chevauchants en un fichier gtf
"""
obFichier1 = open("/home/mbertrand/Bureau/scripts/pipeline/Pic_detection/scripts/uniqueF2Conditions.txt",'r')
output=open("/home/mbertrand/Bureau/scripts/pipeline/Pic_detection/scripts/communEpic.txt","w")

output.write("gene"+"\t"+ "chr" +"\t"+"start"+"\t"+"stop"+"\n")   
for line in obFichier1.readlines():
    ligne=line.split(" ")
      
    chrs=ligne[0]
    chrs=chrs.split("chr")
    chrs=chrs[1]
    start=ligne[1]
    end=ligne[2]
    gene=ligne[3]
    gene=gene.split("\n")
    gene=gene[0]
   
    output.write(gene+"\t"+chrs+"\t"+start+"\t"+end+"\n")
        
