#!/usr/bin/python
# -*- coding: utf-8 -*-

###########
# IMPORTS #
###########
import os.path
import sys
from os import getcwd
import getopt
import re

#########
# Value #
#########

usage = '''

         _______________________________________________________________________________________________________________
        |   Name :                                                                                                      |
        |           compile_counts.py                                                                                   |
        |   But : script qui met en commun les comptages des differentes lib                                            | 
        |                                                                                                               |
        |   EDB (09/07/18)                                                                                              |
        |_______________________________________________________________________________________________________________|
       
'''
   

#############################
# BEGINNING OF THE SOFTWARE #
#############################

## FUNCTIONS

def parse_infile(all_elements_counts, i, infiles):
    """
    all_elements_counts : dico 
    infiles : les fichiers input qui lon passe en arguments
    i : le nombre de infiles

    Parse les differents fichiers pour mettre le comptages en communs

    """
    f1 = file(infiles[i], "r")
   
    for line in f1:
       
            split_line = line.split('\t') #1ligne = 1fenetre
            count = split_line[3] #Le nombre de reads dans chaque fenetre
            element = split_line[0] + "_" + split_line[1] + ":" + split_line[2] 
           
            if (element not in all_elements_counts.keys()):
                all_elements_counts[element] = ["0"]*len(infiles)
            all_elements_counts[element][i] = count
               
   
    f1.close()
    return(all_elements_counts) 
   
   
   
   
def print_data(all_elements_counts, output, infiles):
    """
    Input : output : le nom du fichier de sortie
    infiles : les differents fichiers que lon veut mettre en communs
    all_elements_counts : dico des comptages mis en communs

    Ouput : un fichier au format txt qui contient une entete : element name1 name2 ...
    """
    f2=open(output, "w")
    f2.write("element")
    for f in infiles:
        f2.write("\t" + f.split("_")[0] )
    f2.write("\n")
    for key, values in all_elements_counts.items():
        f2.write(key + "\t" + '\t'.join(values) + "\n")
   
    f2.close()
   
   
   
   

def main(infiles, output):
    """
    Input : infiles sont les fichiers que lon passe en arguments
    output correspond au nom du fichier de sortie (egalement passe en arg dans la ligne de commande)

    Output : Les comptages rassembles pour les differentes lib
    """
    all_elements_counts= {}   
    for i in range (len(infiles)) :
        all_elements_counts = parse_infile(all_elements_counts, i, infiles)
    print_data(all_elements_counts, output, infiles)
    #print_data(genes, go_terms, agrigo)
   
   
   
 
if __name__ == "__main__":

    output_name = "output"

    # getting of parameters
    args = sys.argv[1:]
   
    infiles= []
    # Test for at least 1 argument
    if (args == []):
        print (usage)
        sys.exit()

    try:
        optlist, args = getopt.getopt(args, "ho:", ["help", "output="])

    except getopt.GetoptError, err:
        # print help information and exit:
        print ("option not recognized")
        sys.stderr.write(usage)
        sys.exit(2)

    for opt, arg in optlist:       
        if opt in ("-o", "--output"):
            output = str(arg)       
        elif opt in ("-h","--help") :
            print (usage)
            sys.exit()
        else :           
            assert False, "unhandled option"
            print (usage)
            sys.exit()
   
    infiles = args
    main(infiles, output)

