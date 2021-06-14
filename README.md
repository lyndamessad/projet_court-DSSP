# projet_court-DSSP

# Projet-cours---Assignation-des-structures-secondaires-de-proteines
Réalisation d'un projet court dans le cadre d'une formation universitaire.

	I.	But du script: 
Ce projet consite à implémenter le programme de Désignation de structures secondaires des protéine(DSSP)
Pour le bon fonctionnement des scripts, il faudrait etre dans un environnement conda avec python3.



 II.	Description du programme:
Le programme se base sur le langage python3 sous forme de différentes fonctions.
Il se compose d'un script python exécutable en une ligne de commande sur des fichiers PDB:
    --> python main.py file.pdb
    avec file.pdb pouvant être les fichiers pdb présents dans le dossier Data de ce Github
    Attention : vérifiez que vous ayez bien les fichiers pdb dans le même dossier que les scripts
    --> dans le script main.py, celui importe en module le code pdbHydrogene.py permettant d'exéctuer Reduce sur le fichier PDB entrée dans la commande principale
    



	III.	Les modules et environnement requis:
   • os: le but principal du module OS est d'interagir avec votre système d'exploitation. Très utile pour créer, supprimer ou déplacer des dossiers et parfois de changer ou trouver le répertoire de travail.
   • Sys: le module sys contient des fonctions et des variables spécifiques à l'interpréteur Python lui-même. Ce module est particulièrement intéressant pour récupérer les arguments passés à un script Python lorsque celui-ci est appelé en ligne de commande.
   • Math: le module math est utile pour les manipulations des fonctions et constantes mathématiques de base (sin, cos, exp, pi...)
   • Conda environnement avec python 3 et le programme reduce.3.24.130724 (reduce: version 3.24 07/24/2013, Copyright 1997-2013, J. Michael Word) installé.
   • Getpass : "donne un moyen sécurisé de gérer les invites de mot de passe où les programmes interagissent avec les utilisateurs via le terminal"
   • Numpy : permet d’effectuer des calculs numériques avec Python
   • Argparse : "facilite l'écriture d'interfaces de ligne de commande conviviales, définit les arguments dont il a besoin, et argparse trouvera comment les analyser en dehors de sys. argv
   
   

