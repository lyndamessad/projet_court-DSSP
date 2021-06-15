__auteurs__ = ("Laetitia HOLLIER, Lynda MESSAD DIT MAHTAL")
__version__ = "1.0.0"
__date__ = "2021-05-29"


# Importing the modules necessary for the proper functioning of the script
import sys
import os
import getpass
import argparse
from pdbHydrogene import pdbReduce
import math
from numpy import arccos,sqrt,pi



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
    Return the PDB file given in command
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



def get_trouve_CN(nom):
    """
    Before starting our analysis on proteins, we must check for each of them the distance
    peptide between these residues. If the entire distance found is included
    between 2.2Angstrom and 3Angstrom, then we can continue. Otherwise it will want it
    signify a rupture of the peptide bond and therefore of the carbon skeleton.
    This function takes a pdb file
    Collect the coordinates of the atoms C et N involved in the peptide bond
    """
    with open(nom, "r") as filin1:
        list_C = []
        list_N = []
        for line in filin1:
            dict_C = {}
            dict_N = {}
            i = 1
            if line.startswith("ATOM") and line[12:16].strip() == "C" :
                dict_C["resid"] = i
                dict_C["x"] = float(line[30:38])
                dict_C["y"] = float(line[38:46])
                dict_C["z"] = float(line[46:54])
                list_C.append(dict_C)
            elif line.startswith("ATOM") and line[12:16].strip() == "N" :
                dict_N["resid"] = i
                dict_N["x"] = float(line[30:38])
                dict_N["y"] = float(line[38:46])
                dict_N["z"] = float(line[46:54])
                list_N.append(dict_N)
            i += 1

        return list_C, list_N



def get_distance(xA,yA,zA,xB,yB,zB):
    """
    This function takes 6 floats (coordinates)
    It gives an Euclidean distance between two atoms (Angsstroms)
    returns a float
    """
    xA,yA,zA = float(xA),float(yA),float(zA)
    xB,yB,zB = float(xB),float(yB),float(zB)

    dist = math.sqrt( (xB-xA)**2 + (yB-yA)**2 + (zB-zA)**2 )
    return dist



def get_pre_parsing(pdbH):
    """
    This function takes a PDB file
    Search the atoms C,O,N,H,CA and put the informations in a tuple
    Return a dictionary of lists of tuples
    """
    atoms=["C","O","N","H","CA"]
    with open(pdbH, "r") as filin:
        dict_coord = {}
        for line in filin:
            if line.startswith("ATOM") and line[12:16].strip() in atoms and (line[21:23].strip()=="A" or line[21:23].strip()=="") :
                posi_res = line[22:26].strip()
                atom_name = line[12:16].strip()
                residu_name = line[17:20].strip()
                x = line[31:38]
                y = line[40:46]
                z = line[48:54]

                if posi_res in dict_coord:
                    dict_coord[posi_res].append((residu_name,atom_name,x,y,z))
                else:
                    dict_coord[posi_res]=[(residu_name,atom_name,x,y,z)]

    return dict_coord



def get_parsing(dict_coord):
    """
    This function takes a dictionary of lists of tuples and a PDB file modified by Reduce
    Determinate if each key contains in their list 5 tuples
    return a dico of list of tuples
    """
    dico = {}
    for key in dict_coord:
        if len(dict_coord[key]) >=5 :
            dico[key] = dict_coord[key]
    return dico



def get_energie(r_ON,r_CH,r_OH,r_CN):
    """
    This function takes 4 floats
    It calculates the energie (kcal/mol) of Hbond between two residues (Coulomb's formula)
    returns a float
    # F: dimensional factor
    # Delta : product bteween the unit electron charge (q1 = 0.42e and q2 = 0.20e)
    """
    F = 332
    Delta = 0.084
    E = Delta * ( r_ON**-1 + r_CH**-1 - r_OH**-1 - r_CN**-1 ) * F
    return E




def found_LH(dict_coord):
    """
    This function takes a dictionary of lists of tuples
    Determinate the Hbond, calculate distance and energie
    Return a dictionary with each key is a tuple and the value of the key is a float
    """

    LH_dict ={}
    dic_keys= list(dict_coord.keys())
    dic_values=list(dict_coord.values())
    for i in range(len(dic_keys)-1):
        a = dict_coord[dic_keys[i]] #residu1 
        C = a[1][2:] 
        O = a[2][2:]
        for j in range(len(dic_keys)):
            b = dict_coord[dic_keys[j]]
            N = b[0][2:]
            if len(b) == 5:
                H = b[3][2:]
                r_ON = get_distance(O[0],O[1],O[2],N[0],N[1],N[2])
                r_CH = get_distance(C[0],C[1],C[2],H[0],H[1],H[2])
                r_OH = get_distance(O[0],O[1],O[2],H[0],H[1],H[2])
                r_CN = get_distance(C[0],C[1],C[2],N[0],N[1],N[2])
                energy = get_energie(r_OH,r_CH,r_OH,r_CN)
                if energy < -0.5: 
                    LH_dict[(dic_keys[i],dic_keys[j])] = energy
    return LH_dict



def found_patterns(dict_LH):
    """
    This function takes a dictionary with each key is a tuple and the value of the key is a float
    Determinate the pattern: 3turn, 4turn, 5turn
    Return a dictionary of lists of tuples
    """
    patterns = {} 
    list_turn3 = []
    list_turn4 = []
    list_turn5 = []
    list_none =[]
    list_strand =[]

    for elm in dict_LH:
        a=int(elm[0]) #num de residu
        b=int(elm[1])
        if abs(a-b)==3: #3turns
            list_turn3.append((a,b))
            pat = "3turns"
        elif abs(a-b)==4: #4turns
            list_turn4.append((a,b))
            pat = "4turns"
        elif abs(a-b)==5: #5turns
            list_turn5.append((a,b))
            pat = "5turns"
        elif abs(a-b)<3: #no structure
            list_none.append((a,b))
            pat = "None"
        elif abs(a-b)>5: #brin=strand
            list_strand.append((a,b))
            pat="strand"

        if pat in patterns:
            patterns[pat].append((a,b))
        else:
            patterns[pat]=[(a,b)]
    return patterns



    return helix



def get_bend(dico):
    """
    Takes a first list (dico) having lists (representing each residu) containing dictionaries ((having each atoms of a residu))
    Determine a bend of a residu
    Return a list of integers
    """
    i = 2
    list_bend = []
    dic_keys = list(dico.keys())
    dic_values =list(dico.values())
    while i < (len(dico) -2):
        # Coordinates of CA(i-2)
        x_b,y_b, z_b  = dico[dic_keys[i-2]][0][2], dico[dic_keys[i-2]][0][3], dico[dic_keys[i-2]][0][4]
        # Coordinates of CA(i) 
        x, y, z = dico[dic_keys[i]][0][2], dico[dic_keys[i]][0][3], dico[dic_keys[i]][0][4]
        # Coordinates of CA(i+2)
        x_a, y_a, z_a = dico[dic_keys[i+2]][0][2], dico[dic_keys[i+2]][0][3], dico[dic_keys[i+2]][0][4]
        # Determination of vectors between CA(i-2)-CA(i) and CA(i)-CA(i+2)
        vecteur_CAb_CA = get_vector(x_b,y_b,z_b,x,y,z)
        vecteur_CA_CAa = get_vector(x,y,z,x_a,y_a,z_a)
        # Calculate of the angle between the 2 vectors
        angle = get_angle_bend(vecteur_CAb_CA, vecteur_CA_CAa)
        if angle > 70 : 
            list_bend.append(dic_keys[i])
        i += 1
    return list_bend



def get_vector(x1,y1,z1,x2,y2,z2):
    """
    This function  takes two list of dictionary which contain coordinates
    returns a vector calculating betwwen two points from those two list in a dictionary (floats)
    """
    vecteur = []
    u = float(x1) - float(x2)
    v = float(y1) - float(y2)
    w = float(z1) - float(z2)
    vecteur.append([u,v,w])
    return vecteur



# Helped by https://www.mathweb.fr/euclide/2020/05/15/calculer-la-valeur-dun-angle-avec-le-produit-scalaire/
def get_angle_bend(u,v):
    """
    This function takes 2 dictionary of vectors (float)
    Calculate the angle of a bend
    Return a float
    """
    scalar_product = ( float(u[0][0]) * float(v[0][0]) ) + ( float(u[0][1]) * float(v[0][1]) ) + ( float(u[0][2]) * float(v[0][2]) ) 
    normeU = math.sqrt( float(u[0][0])**2 + float(u[0][1])**2 + float(u[0][2])**2 )
    normeV = math.sqrt( float(v[0][0])**2 + float(v[0][1])**2 + float(v[0][2])**2)
    angle = round(arccos( scalar_product / (normeU * normeV) ) * 180 / pi , 2)
    return angle



# Helped by https://www.mathweb.fr/euclide/2020/05/15/calculer-la-valeur-dun-angle-avec-le-produit-scalaire/
def get_angle_chiral(u,v,w):
    """
    This function takes 3 dictionary of vectors (float)
    Calculate the angle of chirality
    Return a float
    """
    scalar_product = ( float(u[0][0]) * float(v[0][0]) * float(w[0][0]) ) + ( float(u[0][1]) * float(v[0][1]) * float(w[0][1])) + ( float(u[0][2]) * float(v[0][2]) * float(v[0][2])) 
    normeU = math.sqrt( float(u[0][0])**2 + float(u[0][1])**2 + float(u[0][2])**2 )
    normeV = math.sqrt( float(v[0][0])**2 + float(v[0][1])**2 + float(v[0][2])**2 )
    normeW = math.sqrt( float(w[0][0])**2 + float(w[0][1])**2 + float(w[0][2])**2 )
    angle = arccos( scalar_product / (normeU * normeV * normeW) ) * 180 / pi 
    return angle



def get_chirality(dico):
    """
    Takes a first list (list_parsing) having lists (representing each residu) containing dictionaries ((having each atoms of a residu))
    Determine the chirality  
    Return a list of tuple (integers, float)
    """
    i = 1
    list_chirality = []
    dic_keys = list(dico.keys())
    dic_values =list(dico.values())
    while i < (len(dico) -2):
        # Coordinates of CA(i-2)
        x_b,y_b, z_b  = dico[dic_keys[i-2]][0][2], dico[dic_keys[i-2]][0][3], dico[dic_keys[i-2]][0][4]
        # Coordinates of CA(i-) 
        x, y, z = dico[dic_keys[i]][0][2], dico[dic_keys[i]][0][3], dico[dic_keys[i]][0][4]
        # Coordinates of CA(i+1)
        x_a, y_a, z_a = dico[dic_keys[i+1]][0][2], dico[dic_keys[i+1]][0][3], dico[dic_keys[i+1]][0][4]
        # Coordinates of CA(i+2)
        x_af, y_af, z_af = dico[dic_keys[i+2]][0][2], dico[dic_keys[i+2]][0][3], dico[dic_keys[i+2]][0][4]
        # Determination of vectors between CA(i-2)-CA(i) and CA(i)-CA(i+2)
        vecteur_CAb_CA = get_vector(x_b,y_b,z_b,x,y,z)
        vecteur_CA_CAa = get_vector(x,y,z,x_a,y_a,z_a)
        vecteur_CAa_CAaf = get_vector(x_a,y_a,z_a,x_af,y_af,z_af)
        # Calculate of the angle between the 3 vectors
        angle = get_angle_chiral(vecteur_CAb_CA, vecteur_CA_CAa, vecteur_CAa_CAaf)
        if (angle > 0 and angle < 180) or (angle < 0 and angle > -180): 
            list_chirality.append((dic_keys[i], int(angle)))
        i += 1
    return list_chirality




def verify_consec_hel(atoms): 
    """
    This function takes a list of couple of residues containing a tuple
    return those which are consecutive.
    """
    succ=[]
    two=[]

    for i in range(1,len(atoms)-1):
        if int(atoms[i-1][0]+1)==int(atoms[i][0]) :
                    #donc ils se suivent
            if int(atoms[i-1][1]+1)==int(atoms[i][1]) :
                #donc ils sont successifs 
                two.append(atoms[i-1])
                two.append(atoms[i])
                succ.append(two)
                two=[]  
    return succ



def verify_consec_strand(atoms):
    """
    This function takes a list of couple of residues containing a tuple
    return those which are consecutive.
    """
    succ=[]
    two=[]

    for i in range(1,len(atoms)-1):
        if int(atoms[i-1][0]+2)==int(atoms[i][0]) :
                    #donc ils se suivent
            if int(atoms[i-1][1]+2)==int(atoms[i][1]) :
                #donc ils sont successifs 
                two.append(atoms[i-1])
                two.append(atoms[i])
                succ.append(two)
                two=[]  
    return succ


def found_helix(patterns):
    """ 
    This function takes a dictionary of lists of tuples
    If 2 turns are consecutive at position i and i-1 => HELIX
    Return a dictionary of lists of lists of tuples
    """
    dic_helix={}
    helix4 = []
    helix3 = []
    helix5 = []

    two=[]
    succ=[]

    turns=list(patterns.keys())
    for turn in turns:
        if turn == "4turns":
            atoms = list(patterns[turn])
            typ="4-helix"
            helix4 = verify_consec_hel(atoms)
            dic_helix[typ]=helix4
        if turn == "3turns":
            atoms = list(patterns[turn])
            typ="3-helix"
            helix3 = verify_consec_hel(atoms)
            dic_helix[typ]=helix3
        if turn == "5turns":
            atoms = list(patterns[turn])
            typ="5-helix"
            helix5 = verify_consec_hel(atoms)
            dic_helix[typ]=helix5
        if turn=="None":
            atoms = list(patterns[turn])
            typ="None"
            None_str = verify_consec_hel(atoms)
            dic_helix[typ]=None_str
        if turn=="strand":
            atoms = list(patterns[turn])
            typ="strand"
            strand_str = verify_consec_strand(atoms)
            dic_helix[typ]=strand_str

    return dic_helix



def get_dihedral(coor1, coor2, coor3, coor4):
    """
    This function takes 4 coordinates (float)
    Calculate of an dihedre angle  phi,psi
    Return the value of an angle
    """
    import numpy as np

    vec12 = float(coor1[0]) - float(coor2[0]), float(coor1[1]) - float(coor2[1]), float(coor1[2]) - float(coor2[2])
    vec23 = float(coor2[0]) - float(coor3[0]), float(coor2[1]) - float(coor3[1]), float(coor2[2]) - float(coor3[2])
    vec34 = float(coor3[0]) - float(coor4[0]), float(coor3[1]) - float(coor4[1]), float(coor3[2]) - float(coor4[2])

    vec_norm1 = np.cross(vec12, vec23)
    vec_norm2 = np.cross(vec23, vec34)
     
    angle_rad = math.acos(np.dot(vec_norm1, vec_norm2) / (np.linalg.norm(vec_norm1) * np.linalg.norm(vec_norm2)))
    angle_deg = angle_rad * 180 / math.pi

    if np.dot(vec12, vec_norm2) < 0:    #Calcul du produit mixte, pour assigner le bon signe à l'angle.
        return angle_deg
    else:
        return -angle_deg


def get_angles(dict_coord):
    """
    This function takes a dictionary of lists of tuples and a PDB file modified by Reduce
    Determinate an angle
    Return the value of an angle
    """
    coors_N=[]
    coors_C=[]
    coors_CA=[]
    angles1={}
    angles2={}
    angles={}

    residues=list(dict_coord.keys())
    values= list(dict_coord.values())

    for res in residues:
        coors_N.append(dict_coord[res][0][2:5])
        coors_C.append(dict_coord[res][2][2:5])
        coors_CA.append(dict_coord[res][1][2:5])
    #print(coors_N)
    for i in range(0,len(coors_N)-1):
        phi=get_dihedral(coors_C[i-1], coors_N[i], coors_CA[i], coors_C[i])
        psi=get_dihedral(coors_N[i], coors_CA[i], coors_C[i], coors_N[i+1])
        angles1[residues[i]]=phi
        angles2[residues[i]]=psi
    angles["Phi"]=angles1
    angles["Psi"]=angles2
    return angles



def found_helix_angl(dict_res,dico_angles,helix):
    """ pour valider les structures secondaires
        helice alpha : Φ = -60° et Ψ = -50°. +-30
        helice 3.10 : Φ = -49° et Ψ = -26°. +-30
        helice pi : Φ = -57° et Ψ = -70°. +-30
        feuillet β : Φ = -139° et Ψ = +135°. +- 30
        -----
        dict_res : dictionnaires des patern,
        dico_angles des angles ,
        typ : 3-helix/4-helix/5-helix
    """
    helice ={}
    residu =[]
    if str(helix) =='4-helix': #helice alpha  --> H
        for i in range(len(dict_res[helix])):
            for j in dict_res[helix][i]:
                res1=j[0]
                res2=j[1]
                
                psi_1=dico_angles["Psi"][str(res1)]
                phi_1=dico_angles["Phi"][str(res1)]
                
                psi_2=dico_angles["Psi"][str(res2)]
                phi_2=dico_angles["Phi"][str(res2)]
                if (phi_1 < -20 and phi_1 >= -90) and (psi_1 < -20 and psi_1 > -80):
                    residu=[phi_1,psi_1,"H"]
                    helice[res1]=residu
                    residu =[]
                if (phi_2 < -20 and phi_2 >= -90) and (psi_2 < -20 and psi_2 > -80):
                    residu=[phi_2,psi_2,"H"]
                    helice[res2]=residu
                    residu =[]
    if str(helix) =='3-helix' :#helice 3.10helix ---> G
        for i in range(len(dict_res[helix])):
            for j in dict_res[helix][i]:
                res1=j[0]
                res2=j[1]
                
                psi_1=(dico_angles["Psi"][str(res1)])
                phi_1=(dico_angles["Phi"][str(res1)])
                
                psi_2=(dico_angles["Psi"][str(res2)])
                phi_2=(dico_angles["Phi"][str(res2)])
                if (phi_1 < -19 and phi_1 >= -79) and (psi_1 < 4 and psi_1 > -56):
                    residu=[phi_1,psi_1,"G"]
                    helice[res1]=residu
                    residu =[]
                if (phi_2 < -19 and phi_2 >= -79) and (psi_2 < 4 and psi_2 > -56):
                    residu=[phi_2,psi_2,"G"]
                    helice[res2]=residu
                    residu =[]
    if str(helix) =='5-helix': #pi.helix ---> I
        for i in range(len(dict_res[helix])):
            for j in dict_res[helix][i]:
                res1=j[0]
                res2=j[1]
                
                psi_1=(dico_angles["Psi"][str(res1)])
                phi_1=(dico_angles["Phi"][str(res1)])
                
                psi_2=(dico_angles["Psi"][str(res2)])
                phi_2=(dico_angles["Phi"][str(res2)])
                if (phi_1 < -27 and phi_1 >= -87) and (psi_1 < -40 and psi_1 > -100):
                    residu=[phi_1,psi_1,"I"]
                    helice[res1]=residu
                    residu =[]
                if (phi_2 < -27 and phi_2 >= -87) and (psi_2 < -40 and psi_2 > -100):
                    residu=[phi_2,psi_2,"I"]
                    helice[res2]=residu
                    residu =[]
    return helice


def found_strand_angl(dict_res,dico_angles):
    """ 
    to validate secondary structures
    euillet β : Φ = -139° et Ψ = +135°. +- 30
    ----
    dict_res : dictionnaires des patern,
    dico_angles des angles ,
    """
    #strands---> E
    strands ={}
    residu =[]
    for i in range(len(dict_res["strand"])):
        for j in dict_res["strand"][i]:
            res1=j[0]
            res2=j[1]

            psi_1=dico_angles["Psi"][str(res1)]
            phi_1=dico_angles["Phi"][str(res1)]

            psi_2=dico_angles["Psi"][str(res2)]
            phi_2=dico_angles["Phi"][str(res2)]
            if (phi_1 < -109 and phi_1 >= -169) and (psi_1 < 165 and psi_1 > 105):
                residu=[phi_1,psi_1,"E"]
                strands[res1]=residu
                residu =[]
            if (phi_2 < -109 and phi_2 >= -169) and (psi_2 < 165 and psi_2 > 105):
                residu=[phi_2,psi_2,"E"]
                strands[res2]=residu
                residu =[]
    return strands



def get_write(pdbH, dico, bend, chirality, H, G, I, E):
    """
    This function takes the name of a pdb file, a dictionary
    """
    args = str(pdbH).split("H")[0] +"_parsing_results.txt"
    with open(args, "w") as filout:
        filout.write("FILE PARSING OF THE PEPTIDE : {}".format(str(args).split("_")[0]))
        filout.write("\nResidues with coordinates of each of their atoms\n")
        for key in dico:
            filout.write("\n{} : {}\n{} : {}\n{} : {}\n{} : {}\n{} : {}\n".format(key, dico[key][0], key, dico[key][1], key, dico[key][2], key, dico[key][3], key, dico[key][4]))
        filout.write("\n Bend: list of residu\n")
        filout.write("{} \n".format(str(bend)))
        filout.write("\n Chirality: list of residu\n")
        filout.write("{} \n".format(str(chirality)))
        filout.write("\n H (4-helix): list of residu\n")
        filout.write("{} \n".format(str(H)))
        filout.write("\n G (3-helix): list of residu\n")
        filout.write(" {} \n".format(str(G)))
        filout.write("\n I (5-helix): list of residu\n")
        filout.write("{} \n".format(str(I)))
        filout.write("\n E (strand): list of residu\n")
        filout.write("{} \n".format(str(E)))
    return args


def main():
    #arguments=get_args()
    #print(arguments)

    
    print("\n-----------------------------------------------")
    print("              STEPS FOR ARGUMENTS           ")
    print("------------------------------------------------")
    # Verify the existence of arguments given in the python command
    args = verify_args(sys.argv) 
    print(type(args))


    print("")
    print("\n----------------------------------------------")
    print("    STEPS FOR VERIFICATION OF LENGTH PEPTIDE    ")
    print("-----------------------------------------------")
    list_CN = get_trouve_CN(args)
    list_C, list_N, count = list_CN[0], list_CN[1], 0
    distance_liaison_C_N = []
    for i in range(len(list_C) -1):
        distance = get_distance(list_C[i]["x"],list_C[i]["y"], list_C[i]["z"], list_N[i+1]["x"], list_N[i+1]["y"], list_N[i+1]["z"])
        distance_liaison_C_N.append(distance)
        if distance < 3 :
            count += 1
    if len(distance_liaison_C_N) == count:
        print("\n=> The {} peptide bonds are not broken. We can continue".format(count))
    else:
        sys.exit("\n=> ERROR ! Cleavage of peptide bonds !")


    print("")
    print("\n----------------------------------------------")
    print("             STEPS FOR REDUCE   ")
    print("-----------------------------------------------")
    # ---------------------------------------
    # STEPS FOR REDUCE 
    # ---------------------------------------
    # Produce the pdb file which contain the hydrogens atoms added by the Reduce program
    pdbH = pdbReduce(args)
    print("")
    print(pdbH)

    print("")
    print("\n----------------------------------------------")
    print("     STEPS FOR PARSING THE PDB HYDROGEN FILE   ")
    print("----------------------------------------------")
    # Name the file pdbH pre-parsing
    dict_coord = get_pre_parsing(pdbH)
    print("\n        => Step of pre-parsing there were {} résidus".format(len(dict_coord)))
    # Parsing
    dict_coord_parsing = get_parsing(dict_coord)
    print("\n        => After parsing there are {} résidus".format(len(dict_coord_parsing)))


    print("")
    print("\n----------------------------------------------")
    print("        STEPS FOR IDENTIFYING THE PATTERN  ")
    print("----------------------------------------------")
    # Found the Hbond
    LH = found_LH(dict_coord)
    #print(LH)
    patterns = found_patterns(LH)
    #print(patterns)
    
    print("\nRESEARCH OF PATTERNS")
    print("\n=> Found the turn <=")
    # Search of the turns
    # 3turn
    print("\n         => There are {} 3turns".format(len(patterns["3turns"])))
    print("3turns : {}".format(patterns["3turns"]))
    # 4turn
    print("\n         => There are {} 4turns".format(len(patterns["4turns"])))
    print("4turns : {}".format(patterns["4turns"]))
    # 5turn
    print("\n         => There are {} 5turns".format(len(patterns["5turns"])))
    print("5turns : {}".format(patterns["5turns"]))


    print("")
    print("\n----------------------------------------------")
    print("      STEPS FOR IDENTIFYING THE STRUCTURES  ")
    print("----------------------------------------------")
    print("\nRESEARCH OF STRUCTURES")
    
    # characterization of the bend
    bend = get_bend(dict_coord_parsing)
    print("\n         => There are {} CA bend".format(len(bend)))
    print("bend: {}".format(bend))

    # characterization of the chirality
    chirality = get_chirality(dict_coord_parsing)
    print("\n         => There are {} CA chiral".format(len(chirality)))
    print("chirality : {}".format(chirality))
    print("")

    # characterization of the helix
    #get angles phi and psi 
    dico_angles=get_angles(dict_coord)

    # helix from patterns
    dico_paterns=found_helix(patterns)

    # get angles
    dico_angles=get_angles(dict_coord)
    # ---------------- HELIX -------------------
    # 3-helix
    G=found_helix_angl(dico_paterns,dico_angles,'3-helix')

    # 4-helix
    H=found_helix_angl(dico_paterns,dico_angles,'4-helix')

    # 5-helix
    I=found_helix_angl(dico_paterns,dico_angles,'5-helix')

    
    # 4-helix
    print("\n         => There are {} H (4-helix)".format(len(H)))
    print("H : {}".format(H))
    # 3-helix
    print("\n         => There are {} G (3-helix)".format(len(G)))
    print("G : {}".format(G))
    # 5-helix
    print("\n         => There are {} I (5-helix)".format(len(I)))
    print("G : {}".format(I))

    # ---------------- HELIX -------------------
    #found strand

    E= dico_paterns['strand']
    print("\n         => There are {} E (strand)".format(len(E)))
    print("E : {}".format(E))


    print("")
    print("\n----------------------------------------------")
    print("      STEPS FOR WRITINF THE OUTPUT ")
    print("----------------------------------------------")
    write = get_write(pdbH, dict_coord_parsing, bend, chirality, H, G, I, E)


    # ---------------------------------------

if __name__ == '__main__':
    main()