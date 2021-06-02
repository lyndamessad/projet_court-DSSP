__auteurs__ = ("Laetitia HOLLIER, Lynda MESSAD DIT MAHTAL")
__version__ = "1.0.0"
__date__ = "2021-05-29"


# Importing the modules necessary for the proper functioning of the script
import sys
import os
import getpass
import argparse
from pdbHydrogene import pdbReduce

#### script principal à develloper
"""    def getArgs(args):
    
    script prend un argument = fichier pdb  .
    option -r : zapper reduce
    
    file = sys.argv[1:]

     Récupérer les fichiers donnés en arguments à la fonction    if len(args) != 2:
        sys.exit("Nombre d'argument incorrect.\n Le programme prend un fichier pdb. \n")
    elif  file.split(".") != "pdb":
        sys.exit("format de fichier incorrect.\n Le programme ne prend que un fichier type .pdb") 
   
    else:
        if os.path.exists(file):
            print("les fichiers choisis sont: {}\n".format(file))
        else:
            print("le fichier {} n existe pas dans le repertoire{}".format(arg,os.getcwd()))
    return file
"""

def verify_args(args):
    """
    This function takes the arguments given in the python command
    returns the name of the pdb file
    """
    if len(args) != 2:
        sys.exit("Wrong number of arguments.\n"
                 "Usage:  python3 main.py file.pdb")
    f_in = args[1]
    
    extension = f_in[-4:]   #Récupère l'extension du fichier en entrée
    if extension != ".pdb":
        sys.exit("{} does not fit to the .pdb format.\n"
                 "Usage:  python3 main.py file.pdb".format(f_in))
        
    if not os.path.exists(f_in):
        sys.exit("{} does not exist.\n"
                 "Usage: python3 main.py file.pdb".format(f_in))
    return f_in


def get_args():
    """
    COM A FAIRE
    """
    # help
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    # pdb file
    parser.add_argument('-i', '-pdb_file', dest='pdb_file', required=True,
                        help=" Required one pdb file")
    # skype Reduce
    parser.add_argument('-r', '-reduce', dest='reduce', required=False,
                        default=True, action='store_true',
                        help="Add -r to add hydrogene to pdb file")
    parser.add_argument('-no-r', '-no-reduce', dest='no-reduce', required=False,
                        default=False, action='store_true',
                        help="Add -no-r if hydrogene are already added to pdb file")

    return parser.parse_args()


def get_pre_parsing(pdbH):
    """
    This function takes a pdb file which contain the hydrogen atoms given by Reduce.
    The coordinates of the atoms involved in the peptide bond and carbon skeleton are only taken.
    It is written in a pdb file.
    Returns a list of dictionary which have these coordinates and the name of pdb file.
    """
    pdbH_p = str(pdbH).split(".")[0] + "_preparsing." + str(pdbH).split(".")[1]
    with open(pdbH, "r") as filin, open(pdbH_p, "w") as filout:
        list_coord = []
        i = 0 
        line = filin.readlines()
        while i < len(line):
            dict_coord = {}
            if line[i].startswith("ATOM") and (line[i][12:16].strip() == "N" or line[i][12:16].strip() == "C" or line[i][12:16].strip() == "O" or line[i][12:16].strip() == "H" or line[i][12:16].strip() == "CA") :
                filout.write(line[i])
                dict_coord["num_atom"] = str(line[i][6:11].strip())
                dict_coord["atom_name"] = str(line[i][12:16].strip())
                dict_coord["residu_name"] = str(line[i][17:20].strip())
                dict_coord["chain"] = str(line[i][21:22].strip())
                dict_coord["residu"] = int(line[i][22:26].strip())
                dict_coord["x"] = float(line[i][30:38].strip())
                dict_coord["y"] = float(line[i][38:46].strip())
                dict_coord["z"] = float(line[i][46:54].strip())
                list_coord.append(dict_coord)
            elif (line[i].startswith("HETATM") or line[i].startswith("ATOM")) and (line[i][12:16].strip() != "N" or line[i][12:16].strip() != "C" or line[i][12:16].strip() != "O" or line[i][12:16].strip() != "H" or line[i][12:16].strip() != "CA") :
                pass
            else:
                filout.write(line[i])
            i += 1
        return pdbH_p, list_coord


def get_count_residu(liste):
    """
    This function takes a list of dictionary
    Count the number of times a residue is present in the list
    return a dictionary (str,int)
    """
    i = 0
    dico_count = {}
    for i in range(len(liste)) :
        if (str(liste[i]["residu"])) not in dico_count :
            dico_count[str(liste[i]["residu"])] = 1
        else :
            dico_count[str(liste[i]["residu"])] += 1
        i += 1
    return dico_count


def get_deleted(count_residu):
    """
    This function takes a dictionary
    A residu must have at least the 5 atoms (N,CA,C,O,H)
    return a list of integer
    """
    list_del = []
    for NB in count_residu:
        if count_residu[NB] < 5:
            if NB not in list_del :
                list_del.append(NB)
    return list_del


def get_parsing(pre_parsing, residu_del):
    """
    This functions takes 2 arguments
    => pre_parsing : a pdb file 
    => residu_del : list of integer
    Do the parsing
    return the name of a pdb file after the parsing
    """
    parsing = str(pre_parsing).split("_")[0] + "_parsing." + str(pre_parsing).split(".")[1]
    with open(pre_parsing, "r") as filin, open(parsing, "w") as filout:
        i = 0 
        line = filin.readlines()
        while i < len(line):
            if (line[i][22:28].strip() not in residu_del):
                filout.write(line[i])
            i +=1
        return parsing


def get_list(list_pre_parsing, residu_del):
    """
    This function takes 2 arguments:
    => list_pre_parsing : list of dictionary
    => residu_del : list of integer
    Give the coordinates (informations) of the atoms N,CA,C,O,H
    return a list of dictionary
    """
    list_parsing = []
    for i in range(len(list_pre_parsing)):
        if str(list_pre_parsing[i]["residu"]) in residu_del:
            pass
        else:
            list_parsing.append(list_pre_parsing[i])
        i += 1
    return list_parsing


def main():
    #arguments=get_args()
    #print(arguments)
    # ---------------------------------------
    # STEPS FOR ARGUMENTS 
    # ---------------------------------------
    # Verify the existence of arguments given in the python command
    args = verify_args(sys.argv)
    print(type(args))
    # ---------------------------------------
    # STEPS FOR REDUCE 
    # ---------------------------------------
    # Produce the pdb file which contain the hydrogens atoms added by the Reduce program
    pdbH=pdbReduce(args)
    print(pdbH)
    # ---------------------------------------
    # STEPS FOR PARSING THE PDB HYDROGEN FILE 
    # ---------------------------------------
    # Name the file pdbH pre-parsing
    pdbH_pre = get_pre_parsing(pdbH)
    pre_parsing = pdbH_pre[0]
    print("\npre_parsing [0] : {}".format(pre_parsing))
    # Give the pre-parsing list with the coordinates of the atoms N,CA,C,O,H
    list_pre_parsing = pdbH_pre[1]
    print("\nlist_pre_parsing [0] : {}".format(list_pre_parsing[0]))
    
    # Count of atoms per residu 
    count_residu = get_count_residu(list_pre_parsing)
    #print("\ncount_residu : {}".format(count_residu))
    # If for a residu, there are less than 5 atoms, this residu must be deleted because it means that the 5 atoms (N,CA,C,O,H) aren't there !
    residu_del = get_deleted(count_residu)
    print("\n The residus to delete (residu_del) : {}".format(residu_del))

    # Parsing of the protein 
    pdbH_parsing = get_parsing(pre_parsing, residu_del)
    print("\npdbH_parsing : {}".format(pdbH_parsing))
    # Give the parsing list with the coordinates of residu contaning all of these atoms: N,CA,C,O,H
    list_parsing = get_list(list_pre_parsing, residu_del)
    print("\nlist_parsing [0] : {}".format(list_parsing[0]))
     # ---------------------------------------

if __name__ == '__main__':
    main()
