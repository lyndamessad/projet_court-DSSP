"""
le script pdbHydrogene est conçu pour des fichiers au format .pdb. 
En entrant le nom d'un fichier pdb d'une protéine déja existant dans le
repertoire, le module récupere le fichier associé, excécute le programme Réduce 
(préalablement installé sur la machine) afin d'ajouter les liaisons hydrogenes. 
Un fichier output est founit pour chaque protéine en input 
"""
import sys
import os


def pdbReduce(pdb):
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