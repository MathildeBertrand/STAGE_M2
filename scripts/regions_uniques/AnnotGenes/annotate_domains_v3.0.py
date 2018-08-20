#!/usr/bin/python
# -*- coding: utf-8 -*-

###########
# IMPORTS #
###########
import os.path
import sys
from os import getcwd
import getopt

#########
# Value #
#########

usage = '''

         _______________________________________________________________________________________________________________
        |   Name :                                                                                                      |
        |           annotate_domains_v3.0.py                                                                            |
        |                                                                                                               |
        |   Synopsis :                                                                                                  |
        |           python annotate_domains_v3.0.py -f <INPUT> -a <annot> -e <extension> -t <THRESHOLD> -o <OUTPUT> -1  |
        |                                                                                                               |
        |   Description :                                                                                               |
        |          annotate_domains_v3.0.py script permits to annotate some domains given in input in gff format or     |
        |       in  text format like : the header IS NEEDED !                                                           |
        |               domain   omosome      Start   Stop                                                            |
        |               dom1     1    7173328 7173574                                                                 |
        |               dom2     4    1339100 140300                                                                  |
        |               dom3     4    1539100 160300                                                                  |
        |                                                                                                               |
        |       This script permits to give your own annotation file. It gives the list of annotations which :          |
        |               * overlaps                                                                                      |
        |               * or are at a distance shorter than a threshold                                                 |
        |                                                                                                               |
        |   DIFFERENCE WITH VERSION 2 : Use your own annotation file                                                    |
        |                                                                                                               |
        |   Options :                                                                                                   |
        |               -f, --infile (HEADER NEEDED - no fixed text for the header)                                     |
        |                       file with ID, start et stop of the needed region (REQUIRE).                             |
        |                       ==> txt file please !                                                                   |
        |                                                                                                               |
        |                       gene          start   stop                                                            |
        |                       AT5G49160        5    19948956        19949456                                        |
        |                       AT5G39550        5    15857866        15858366                                        |
        |                                                                                                               |
        |               -a, --annotFile                                                                                 |
        |                       Annotation file (REQUIRE)                                                               |
        |                       Must be a GFF file with ID for the first attribute in col 9                             |
        |                                                                                                               |
        |               -e, --extension                                                                                 |
        |                       format of your input file - default (txt)                                               |
        |                       allowed values : TXT / GFF                                                              |
        |                       Beware : with GFF format, the name of your domains MUST be in the attribute field       |
        |                               in the first position.                                                          |
        |                                                                                                               |
        |               -t, --threshold                                                                                 |
        |                       the needed threshold (NOT NEEDED) - default = 1000 (1kb)                                |
        |                       A negative threshold represents the minimal distance for overlaps.                      |
        |                                                                                                               |
        |               -1, --1line                                                                                     |
        |                       only 1 line per domains (NOT NEEDED) - defaut 1 lin per domain                          |
        |                                                                                                               |
        |               -o, --output                                                                                    |
        |                       output (NOT NEEDED) - defaut "output.txt"                                               |
        |                                                                                                               |
        |                                                                                                               |
        |   Notations :                                                                                                 |
        |               B - the element is BEFORE the domain                                                            |
        |                        _______      _______                                                                   |
        |                       |  ELT  |    |  DOM  |                                                                  |
        |                       |_______|    |_______|                                                                  |
        |                                                                                                               |
        |               A - the element is AFTER the domain                                                             |
        |                        _______      _______                                                                   |
        |                       |  DOM  |    |  ELT  |                                                                  |
        |                       |_______|    |_______|                                                                  |
        |                                                                                                               |
        |               OA - overlap with the element AFTER the domain                                                  |
        |               OB - overlap with the element BEFORE the domain                                                 |
        |               DIE - domain included into the element                                                          |
        |               EID - element included into the domain                                                          |
        |                                                                                                               |
        |   Exemple: python annotate_domains_v3.0.py -f MY_doms.txt -t 500                                              |
        |                                                                                                               |
        |   ED (05/03/12)    adaptee en juin 2018 (Mathilde Bertrand) pour permettre annotation des pics 
            avec une annotation de genes                                                                                |
         _______________________________________________________________________________________________________________


'''



#############
# FUNCTIONS #
#############

        
def bigger_elt_lg(table):
    # donne la longueur du plus grand element d'un tableau
    val = 0
    for elt in table:
       size = len(elt)
       if (size > val):
           val = size

    return val

    
def get_ (infile, format_infile):
    """
    input : le fichier et le format
    output : une liste des  s pst dans le fichier
    Principe : Le fichier est trie par ordre des  s et ne retient le nom des  s que une et une seule fois dans un fichier intermediaire
    Le fichier est ensuite lu pour stocker le nom des  s dans une liste. Il est ensuite supprime
    """
     omosomes = []
    
    os.system('''awk -F "\t" '{print $2}'  ''' + infile + " | sort | uniq > infile_ .txt")
    os.system("rm infile1.txt")
    
    f1 = open("infile_ .txt", "r") #le fichier trie par ordre de  s
    
    while True:
        line = f1.readline()
        
        if (line == ""):
            break
            
        split_line = line.split("\t")
         omo = split_line[0].strip()
         omosomes.append( omo)
          
    f1.close()
    os.system("rm infile_ .txt")
    return  omosomes 


def get_data( omosomes, infile, format_infile, thres, annot_file):
    """
    input :  omosomes : une liste des  s pst dans linput
    """

    # data about domains
    domain_ttl = {}
    domains_all = {}
    
    for  omo in  omosomes: #Pour chaque  s de la liste
        print("\n\t management of  s number "+ omo)
       
        # management of annotation data
        #On recupere toutes les annotations sur le  s qui nous interesse puis on tri le fichier. Enfin, on supprime le fichier non trie
        os.system("grep -w '^"+ omo+"' "+ annot_file + " > annot_"+ omo) 
        os.system("sort -k 4 -n annot_"+ omo + " > annot_S_"+ omo)
        annot_file_sort = 'annot_S_'+ omo #Fichier dannotation intermediaire
        os.system("rm annot_"+ omo)
        
        
        # management of infile
        os.system("grep -w "+ omo+" "+infile+ " > infile_"+ omo)
        if (format_infile == "txt"):
            os.system(''' awk -F "\t" '{print $1 "\t" $3 "\t" $4}' infile_''' +  omo + " > infile_out_"+ omo)
        else:
            os.system(''' awk -F "\t" '{split($9, a, ";"); split(a[1], b, " "); print b[2] "\t" $4 "\t" $5}' infile_''' +  omo + " > infile_out_"+ omo)
            
        os.system("sort -k 2 -n infile_out_"+ omo + " > infile_S_"+ omo)
        region_file = 'infile_S_'+ omo #fichier intermediaire qui contient les infos des pics par regions (les unes a la suite des autres)
        os.system("rm infile_"+ omo)  
        os.system("rm infile_out_"+ omo)         
        
        # DATA comparison
        data_domain = {}
            
        # rescue of infile data in a dictonnary
        #########################################
        f1 = open(region_file, 'r')
        while True:
            line_infile = f1.readline()
            
            if (line_infile == ''):
                break
            
            split_line = line_infile.split("\t")
            
            rg = split_line[0]
            start_rg = (split_line[1])
            stop_rg = (split_line[2].strip())
         
            
            if (start_rg != "start" and stop_rg !="stop" and rg !="gene"):
                start_rg=int(start_rg)
                stop_rg=int(stop_rg)
                
                data_domain[start_rg] = [rg, stop_rg]
                domain_ttl[rg] = [ omo, start_rg, stop_rg]
            
        print ("\t\t getting of region data for "+ omo+" ------ ok") 
            
            
        # comparison with all_annot
        ##########################
        f2 = open(annot_file_sort, 'r')
        
        starts = (data_domain.keys())
        starts.sort()    
        
       
        list_all_annot = {}       
        while True:
            line_all_annot = f2.readline()
      
            if (line_all_annot == ''):
                break

            split_all_annot = line_all_annot.split("\t")
            
           
            start_annot = int(split_all_annot[3])
            stop_annot = int(split_all_annot[4])
            strand_annot = str(split_all_annot[6])
            test_ID = len(split_all_annot[8].strip())
            if (test_ID > 1):
                #ID_annot = (split_all_annot[8].strip()).split(";")[0].split("ID=")[1]
                ID_annot=split_all_annot[8].strip()
              
                
            else:
                #print(split_all_annot[8])
                #ID_annot = (split_all_annot[8].strip()).split("gene_id")[1]
                ID_annot=split_all_annot[8].strip()

           
            cpt_starts = 0
     
            for i in range(cpt_starts, len(starts)): 
    
                start_current_rg = int(starts[i])
                atg_current = data_domain[start_current_rg][0]
                stop_current_rg = int(data_domain[start_current_rg][1])     
       
                
                if (domains_all.get(atg_current)):
                    data_old_annot = domains_all[atg_current]
                 
                else: 
                    data_old_annot = []
                    
                # CAS 1 : complete overlap 
                ### CAS a : the domain is totally included into the element i.e. the annot
                if ((start_annot <= start_current_rg) and (stop_annot >= stop_current_rg)):
                    data = ["annot", ID_annot, start_annot, stop_annot, strand_annot, "DIE", '*'] # DIE for domain in element   
                    data_old_annot.append(data)
                    domains_all[atg_current] = data_old_annot
                    
                ### CAS b : the element is totally included into the domain 
                elif ((start_current_rg <= start_annot) and (stop_current_rg >= stop_annot)):
                    data = ["annot", ID_annot, start_annot, stop_annot, strand_annot, "EID", "*"] # EID for element in domain
                    data_old_annot.append(data)
                    domains_all[atg_current] = data_old_annot
                
                # CAS 2 : non complete overlap 
                elif ( ((start_annot < start_current_rg) and (stop_annot >= start_current_rg) and (stop_annot <= stop_current_rg))
                       or ((start_annot >= start_current_rg) and (start_annot <= stop_current_rg) and (stop_annot > stop_current_rg)) ):
                    if ((start_annot < start_current_rg) and (stop_annot >= start_current_rg) and (stop_annot <= stop_current_rg)):
                        dist = stop_annot - start_current_rg
                        annot = "OB" # OA for Overlap Before
                    else :
                        dist = stop_current_rg - start_annot
                        annot = "OA" # OA for Overlap After

                    if (thres <0):
                        if (dist >= -thres):
                            data = ["annot", ID_annot, start_annot, stop_annot, strand_annot, annot, dist]
                            data_old_annot.append(data)
                            domains_all[atg_current] = data_old_annot
                    else :
                        data = ["annot", ID_annot, start_annot, stop_annot, strand_annot, annot, dist]
                        data_old_annot.append(data)
                        domains_all[atg_current] = data_old_annot
                
                # CAS 3 : le annot est AVANT - dans la limite du seuil
                elif ((stop_annot < start_current_rg) and (start_current_rg - stop_annot < thres )) :
                    dist = start_current_rg - stop_annot 
                    data = ["annot", ID_annot, start_annot, stop_annot, strand_annot, "B", dist] # B for BEFORE domain
                    data_old_annot.append(data)
                    domains_all[atg_current] = data_old_annot
                    
                # CAS 4 : le annot est APRES - dans la limite du seuil
                elif ((start_annot > stop_current_rg) and (start_annot - stop_current_rg < thres )):
                    dist = start_annot - stop_current_rg 
                    data = ["annot", ID_annot, start_annot, stop_annot, strand_annot, "A", dist] # A for AFTER domain   
                    data_old_annot.append(data) 
                    domains_all[atg_current] = data_old_annot
                    
                else :
                    cpt_starts = cpt_starts + 1
                
            
        print ("\t\t getting of annotation data ------ ok") 
                

                
        os.system("rm annot_S_"+ omo)
        os.system("rm infile_S_"+ omo)
       

                    
        
        f1.close()
        f2.close()

    os.system("rm infile.txt") 
    return(domain_ttl, domains_all)
  

# Create outfile
#****************
def write_data(domain_ttl, domains_all, outfile, only_one_line):
    
    if (only_one_line):
        print("\n\t DATA PRINTING")

        values_all = list(domains_all.values())
        val_all = bigger_elt_lg(values_all)
        f1 = open(outfile, "w")
        
        title_all = "DOMAIN\t OMOSOME\tSTART_DOMAIN\tSTOP_DOMAIN"
        for u in range(val_all):
            title_all = title_all + "\tELEMENT_"+ str(u+1) + "_NAME\tELT_type\tELT_start\tELT_stop\tELT_strand\tposition_ELT\tdistance"
        f1.write(title_all+"\n")

        keys = list(domain_ttl.keys())
        keys.sort()
        for key in keys :
            if (domains_all.get(key)):
                all_data = domains_all[key]
                line = key + "\t" + str(domain_ttl[key][0]) + "\t" +  str(domain_ttl[key][1]) + "\t" +  str(domain_ttl[key][2])
                data = ""
                for i in range(len(all_data)):
                    data_allType = all_data[i]
                    
                    for dat in data_allType :
                        data = data + "\t" + str(dat)
                        
                f1.write(line + data + "\n")
            else :
                f1.write(key + "\t" + str(domain_ttl[key][0]) + "\t" +  str(domain_ttl[key][1]) + "\t" +  str(domain_ttl[key][2]) + "\tNO_DATA\n")
                
        f1.close()
        print("\nYour output files are: " + outfile + "\n")


    else:
        print("\n\t DATA PRINTING")

        values_all = list(domains_all.values())
        val_all = bigger_elt_lg(values_all)            
        f2 = open(outfile, "w")

        title_all = "DOMAIN\t OMOSOME\tSTART_DOMAIN\tSTOP_DOMAIN\tEL_TYPE\tELT_NAME\tELT_stop\tELT_start\tELT_strand\tPOSITION\tDISTANCE"
        f2.write(title_all+"\n")

        keys = list(domain_ttl.keys())
        keys.sort()
        for key in keys :
            if (domains_all.get(key)):
                all_data = domains_all[key]
                #line = key + "\t" + str(domain_ttl[key][0]) + "\t" +  str(domain_ttl[key][1]) + "\t" +  str(domain_ttl[key][2])
                for i in range(len(all_data)):
                    data_allType = all_data[i]
                    data = ""
                    for dat in data_allType:
                       
                        data = data + "\t" + str(dat)
                    f2.write(key + "\t" + str(domain_ttl[key][0]) + "\t" +  str(domain_ttl[key][1]) + "\t" +  str(domain_ttl[key][2]) + data + "\n")
                       
                            
            else :
                f2.write(key + "\t" + str(domain_ttl[key][0]) + "\t" +  str(domain_ttl[key][1]) + "\t" +  str(domain_ttl[key][2]) + "\tNO_DATA\n")
                
        f2.close()
        print("\nYour output file is: " + outfile + "\n")
            



#############################
# BEGINNING OF THE SOFTWARE #
#############################

def main(infile, annotFile, format_infile, thres, outfile, only_one_line):

   # good_files = test_files(infile, format_infile)
    
     omosomes = get_ (infile, format_infile)            
    domain_ttl, domains_all = get_data( omosomes, infile, format_infile, thres, annotFile)
    write_data(domain_ttl, domains_all, outfile, only_one_line)
    if (os.path.exists("infile.txt")):
        os.system("rm infile.txt")




if __name__ == "__main__":

    outfile = "outfile.txt"
    
    # getting of parameters
    args = sys.argv[1:]
    thres = 1000
    extension = "txt"
    only_one_line = False

    # Test for at least 1 argument
    if (args == []):
        print (usage)
        sys.exit()

    try:
        optlist, args = getopt.getopt(args, "hf:a:e:t:o:1", ["help", "file=", "annotFile", "extension", "threshold=", "output=", "1line"])

    except getopt.GetoptError, err:
        # print help information and exit:
        print ("option not recognized")
        sys.stderr.write(usage)
        sys.exit(2)

    for opt, arg in optlist:        
        if opt in ("-f", "--file"):
            infile= str(arg)      
        elif opt in ("-a", "--annotFile"):
            annotFile= str(arg)
        elif opt in ("-e", "--extension"):
            extension = str(arg)
        elif opt in ("-t", "--threshold"):
            thres = int(arg)
        elif opt in ("-o", "--output"):
            outfile = str(arg)
        elif opt in ("-1", "--1line"):
            only_one_line = True
        elif opt in ("-h","--help") :
            print (usage)
            sys.exit()
        else :
            assert False, "unhandled option"
            print (usage)
            sys.exit()


    # checking of needed parameters
    if (infile == ""):
        assert False, "no infile specified"
        print(usage)
        sys.exit()

    format_infile = extension.lower()
    if ((format_infile != "txt") and (format_infile != "gff")):
        print(usage)
        print("\nProblem with the infile format : " + str(format_infile) + "\nThe system managed txt and gff files uniquely")
        sys.exit()

    
    


    main(infile, annotFile, format_infile, thres, outfile, only_one_line)

