#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
#Mathilde Bertrand
#7 mai 2018
#Les differentes fonction qui sont appellees dans le script compare_peak.py
"""

def OverlapByTwo(dico):
    """
        Regarde les pics 2 a 2 et definit si il y a des overlap
        input : un dico contenant pour chaque reads : le type = fichier1 ou fichier2, le nom du pic, le numero du  s et les positions sur le  s. Contient aussi la longueur des pics
        
        Output : retourne 3 dicos :
        -overlapByTwo : un dico doverlappe 2 a 2 qui pour chaque pic contient son nom en cle, sa position sur le  s
        -temp : les pics qui sont overlappants
        -deja_vu : lensemble des pics que len recontre
        => Les pics qui sont dans deja_vu mais pas dans temp sont des pics uniques
    """
   
    rep1=0
    rep2=0
    deja_vu={} #Dico qui permet de stocker tous les noms rencontres
    temp={} #Dico qui permet de stocker tous les noms des pics qui sont overlappants

    overlapByTwo={} #Dico des overlappes 2 a 2

    for cle1,valeur in dico[1].items():
        names1=cle1
        chrs1=valeur[0]
        start1=valeur[1]
        end1=valeur[2]
       
    
        for cle2,valeur in dico[2].items():
            names2=cle2
            chrs2=valeur[0]
    
            #on stocke les noms 
            deja_vu[names1]=[] 
            deja_vu[names2]=[]
        
            if (chrs1== chrs2): #Les pics des deux fichiers sont sur le meme  s
                start2=valeur[1]
                end2=valeur[2]
                
                startBis=""
                endBis=""
                
               
                #Test des differentes conditions d overlappe :
                if (end1<end2 and start1>start2): #le pic1 est compris dans le pic2 => le pic2 est donc le plus grand
                    startBis=start2
                    endBis=end2
                elif (end1>end2 and start2>start1): #Le pic2 est compris dans le pic1 => le pic1 est donc le plus grand
                    startBis=start1
                    endBis=end1
                elif (end2>start1 and start2<start1 and end1 >= end2) : #Le pic 2 overlap mais deborde sur la gauche =>on retient donc le start du pic2 et le end du pic1 pour avoir le plus grand overlap possible
                    startBis = start2
                    endBis = end1
                elif (end1>start2 and start1<=start2 and end1 < end2)  :# le pic1 overlapp et deborde sur la gauche => start du pic1 et end du pic2 ou Le pic2 overlapp mais deborde sur la droite => start du pic1 et end du pic2
                    startBis = start1
                    endBis = end2
                                       
              
                #Remplissage dico des overlappes 2 a 2 :    
                if (startBis!="" and endBis!=""): #Si les deux pics se chevauchent
                    
                    temp[names1]=[] #on les stocke dans un dico temporaire (permettra den deduire les pics non chevauchants)
                    temp[names2]=[]
                    names=names1+"-"+names2 #Concanteantion des noms
                    overlapByTwo[names]=[]
                    overlapByTwo[names].append(chrs1)
                    overlapByTwo[names].append(startBis)
                    overlapByTwo[names].append(endBis)
                    overlapByTwo[names].append(int(endBis)-int(startBis))
                
                    rep1=rep1+1
                    rep2=rep2+1
                
    return(overlapByTwo,temp,deja_vu)

def OverlapMutliple(overlapByTwo):
    """
    Pour chaque overlappe 2 a 2 definit plus tot, on regarde si il y a un overlappe multiple : exemple si A overlappe B et B overlappe C
    Input : le dico de chevauchement  2 a 2
    Output : un dico de chevauchement multiple de pics 
    """
    comm={} #Dictionnaire doverlappe multiple

    for cle,valeur in overlapByTwo.items():
        
        name=cle
        chrs=valeur[0]
        start=valeur[1]
        stop=valeur[2]
        
        if  chrs not in comm: #Le  chrs na jamais ete rencontre 
            comm[chrs]={}
        comm[chrs][name]=[]
        comm[chrs][name].append(start)
        comm[chrs][name].append(stop)
        
        for cle,valeur in comm.items():
            for noms,position in valeur.items():
                
                if cle==chrs and noms!=name: #On est sur le meme chrs mais pas sur le meme pic
                    i=0
                    complexe=[]
                    starte=[]
                    stoppe=[]
                    
                    
                    if(int(position[0]) <= int(stop) <= int(position[1]) and int(position[0]) <= int(start) <= int(position[1])): #Le nouveau pic inclut celui qui etait present
                           #On modifie le nom
                           #On supprime les deux pics detectes dans le dico
                           #On cree une nouvelle entree dans le dico qui contient le nouveau nom et les bornes corrigees
                        newName=name+"-"+noms
                                                                           
                        if noms in comm[chrs]:
                            del comm[chrs][noms]
                        if name in comm[chrs]:
                            del comm[chrs][name]
                            
                        comm[chrs][newName]=[]
                        comm[chrs][newName].append(position[0])
                        comm[chrs][newName].append(position[1])
                        i=i+1
                        complexe.append(newName)
                        starte.append(position[0])
                        stoppe.append(position[1])
                
                    elif(int(start) <= int(position[1]) and int(position[1]) <= int(stop) and int(start) <= int(position[0]) and int(position[0])<= int(stop)):#Le nouveau pic est inclut dans celui qui etaient deja present 
                        i=i+1
                      
                        #Alors on change le nom et on fixe les nouvelles bornes
                        newName=name+"-"+noms
                        if noms in comm[chrs]:
                            del comm[chrs][noms]
                        if name in comm[chrs]:
                            del comm[chrs][name]
                            
                        comm[chrs][newName]=[]
                        comm[chrs][newName].append(start)
                        comm[chrs][newName].append(stop)
                        complexe.append(newName)
                        starte.append(start)
                        stoppe.append(stop)
                        
                    elif(int(start) <= int(position[1]) <= int(stop) and int(position[0])<int(start)): #Le nouveau pic est inclut et deborde sur la droite
                        #Alors on change les bornes en prenant le min et le max et on change le nom
                        i=i+1
                        
                        newName=name+"-"+noms
                        if noms in comm[chrs]:
                            del comm[chrs][noms]
                        if name in comm[chrs]:
                            del comm[chrs][name]
                        
                        comm[chrs][newName]=[]
                        comm[chrs][newName].append(position[0])
                        comm[chrs][newName].append(stop)
                            
                        complexe.append(newName)
                        starte.append(position[0])
                        stoppe.append(stop)
                        
                    elif(int(position[0]) <= int(stop) <= int(position[1]) and int(start)<int(position[0])):#On deborde sur la gauche
                        i=i+1
                        
                        newName=name+"-"+noms
                        complexe.append(newName)
                        
                        if noms in comm[chrs]:
                            del comm[chrs][noms]
                        if name in comm[chrs]:
                            del comm[chrs][name]
                          
                        comm[chrs][newName]=[]
                        comm[chrs][newName].append(start)
                        comm[chrs][newName].append(position[1])
                        starte.append(start)
                        stoppe.append(position[1])
                    
                       
    return(comm) #Retourne un dico doverlappe multiple

def writeOutput(outputName,OutputOverlapMutliple,diff,dico):
    """
    Input : 
    - outputName : Le nom qui sera contenu dans tous les fichiers de sortie
    - OutputOverlapMutliple : la sortie de la fonction OverlapMutliple = un dico d'overlappe
    -  diff : un dico des pics uniques
    - dico : le dico de tous les pics
    
    output : les fichiers de sortie au format txt
    """
                
    output=open("Everything"+outputName, "w")
    uniques=open("uniqueF1"+outputName,"w") 
    uniques1=open("uniqueF2"+outputName,"w") 
    commun=open("commun"+outputName,"w") 
    
    comm={}
    comm=OutputOverlapMutliple

    for cle,valeur in comm.items():
        for nom,position in valeur.items():
        
            commun.write("chr"+str(cle)+"\t"+position[0]+"\t"+position[1]+"\t"+nom+"\n")
            output.write("chr"+str(cle)+"\t"+position[0]+"\t"+position[1]+"\t"+nom+"\n")
    
    i=0
    while i < len(diff):
        name=list(diff)[i]
        i=i+1
    
        for cle1,valeur in dico[1].items():
            names1=cle1
            if (names1==name):
                chrs1=valeur[0]
                start1=valeur[1]
                end1=valeur[2]
                l1=valeur[3]
                uniques.write("chr"+chrs1+"\t"+start1+"\t"+end1+"\t"+ names1+"\n")
                output.write("chr"+chrs1 +"\t"+start1+"\t"+end1+"\t"+ names1+"\n")
            
        for cle2,valeur in dico[2].items():
            names2=cle2
            if (names2==name):
                chrs2=valeur[0]
                start2=valeur[1]
                end2=valeur[2]
                l2=valeur[3]
                uniques1.write("chr"+chrs2 +"\t"+start2+"\t"+end2+"\t"+ names2+"\n")
                output.write("chr"+chrs2 +"\t"+start2+"\t"+end2+"\t"+ names2+"\n")

    output.close()
    uniques.close()
    uniques1.close()
    commun.close()

    
