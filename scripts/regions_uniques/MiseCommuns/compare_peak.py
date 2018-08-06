#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Mathilde Bertrand
7 mai 2018

But : determiner pour deux fichiers de pics ceux qui sont communs au deux et ceux qui sont uniques

Input : les deux fichiers de pics EPIC

Output : 
4 fichiers de sortie au format txt avec pour chaque pic  : son nom, son  s, sa position sur le  s
- 1 fichier global : contenant les pics communs aux replicats et ceux qui ne sont pas communs
- 1 fichier contenant les pics uniques a un fichier
- 1 fichier contenant les pics uniques a lautre fichier
- 1 fichier contenant que les pics communs

Ligne de commande : python compare_peak.py -f1 file1 -f2 file2 -o output.txt

"""

import functions as func

if __name__ == "__main__":
    import argparse,os,sys
    
########################################################
#                   Arguments
########################################################
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f1", help="Le premier fichier de pics")
    parser.add_argument("-f2", help="le deuxieme fichier de pics")
    parser.add_argument("-o", help="nom du fichier de sortie")
  
    args = parser.parse_args()

    if args.f1 == None or args.f2 == None or args.o==None : #si l'un des arguments est vide
        parser.print_help()
        parser.exit()
    
    file1=os.path.abspath(args.f1)
    file2=os.path.abspath(args.f2)
    output=os.path.basename(args.o)
 

########################################################
#                   Lecture des fichiers
########################################################

    dic={} 
    i=1
  
    while i<=2:
        typ=i
        dic["type"]=[]
        dic["type"].append(typ) 
        dic[typ]={}
    
        if(i==1):
            obFichier1 = open(file1,'r') #fichier1
        else:
            obFichier1 = open(file2,'r') #fichier2
        
        j=0
        for line in obFichier1.readlines(): #Pour chaque lignes du fichier
            j=j+1 #Le numero de peak
            ligne=line.split(" ")
                    
            chrs=ligne[1]
            chrs=chrs.split('"')
            chrs=chrs[1]
            chrs=chrs.split("chr")
            chrs=chrs[1]
         
          
            start=ligne[2]
            end=ligne[3]
           
                
   
            if (len(ligne)==10):
                name="Peak"+str(i)+"_"+str(j)
               
                
            else:
                name=ligne[3]
                name=name.split("\n")
                name=name[0]
        
            lenght=int(end)-int(start) #la longueur du pic
        
            if (name not in dic):
                dic[typ][name]=[]
                dic[typ][name].append(chrs)
                dic[typ][name].append(start)
                dic[typ][name].append(end)
                dic[typ][name].append(lenght)
        
        obFichier1.close()
        i=i+1

    overlappe=func.OverlapByTwo(dic)
    
    #Recherche des overlappes multiples
    
    multiples = func.OverlapMutliple(overlappe[0])
 
    #Recherche des pics uniques a un fichier :
    difference=set(overlappe[2]) - set(overlappe[1])
   
   
########################################################
#           Ecriture des fichiers de sortie
########################################################

    func.writeOutput(output,multiples,difference,dic)
