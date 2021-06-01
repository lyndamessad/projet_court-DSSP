"""
The pdbHydrogene script is designed for files in .pdb format.
By entering the name of a pdb file of a protein already existing in the
directory, the module retrieves the associated file, runs the Reduce program
(previously installed on the machine) in order to add the hydrogen bonds.
An output file is provided for each protein in input
"""
# Importing the modules necessary for a proper functioning of the script
import sys
import os


def pdbReduce(pdb):
    """
    This function takes a pdb file
    returns a pdb file which has the hydrogen atoms added by the Reduce program
    """
    if os.popen("reduce -version") == "/bin/sh: 1: reduce: not found":
        print(" reduce n est pas installer. \n")
        print("Merci d installer reduce : conda install -c mx reduce")

    else:
        argH = str(pdb).split(".")[0] +"H." + str(pdb).split(".")[1]
        os.system(" reduce {} > {}".format(pdb, argH))
    return argH
    
#programme principal:
if __name__ == "__main__":
    #liste_arg = getArgs(sys.argv)
    #pdbReduce(liste_arg)
    pass
