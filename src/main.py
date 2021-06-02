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


def main():
    #arguments=get_args()
    #print(arguments)
    
    # Verify the existence of arguments given in the python command
    args = verify_args(sys.argv)
    print(type(args))
    # Produce the pdb file which contain the hydrogens atoms added by the Reduce program
    pdbH=pdbReduce(args)
    print(pdbH)
    # ---------------------------------------
    # STEPS FOR PARSING THE PDB HYDROGEN FILE 
    # ---------------------------------------
    # Name of the file pdbH pre-parsing
    pdbH_pre = get_pre_parsing(pdbH)
    pre_parsing = pdbH_pre[0]
    print("\npre_parsing : {}".format(pre_parsing))
    # List pre-parsing with the coordinates of the atoms N,CA,C,O,H
    list_pre_parsing = pdbH_pre[1]
    print("\nlist_pre_parsing : {}".format(list_pre_parsing[0]))
    
    # Count of atoms per residu 
    count_residu = get_count_residu(list_pre_parsing)
    print("\ncount_residu : {}".format(count_residu))
    # If for a residu, there are less than 5 atoms, this residu must be deleted because it means that the 5 atoms (N,CA,C,O,H) aren't there !
    residu_del = get_deleted(count_residu)
    print("\nresidu_del : {}".format(residu_del))

    # Parsing of the protein 
    pdbH_parsing = get_parsing(pre_parsing, residu_del)
    print("\npdbH_parsing : {}".format(pdbH_parsing))
    # List parsing with the coordinates of residu contaning all of these atoms: N,CA,C,O,H
    list_parsing = get_list(list_pre_parsing, residu_del)
    print("\nlist_parsing : {}".format(list_parsing[0]))
     # ---------------------------------------

if __name__ == '__main__':
    main()
