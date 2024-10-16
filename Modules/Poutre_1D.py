# Unites : N-mm-MPa

import numpy as np
import math

from Modules.Fonctions_partagées import (calculer_k_poutre1d, assembler_matrice, extraire_matrice, extraire_vecteur,
                                         reconstruire_vecteur, assembler_vecteur)

# Variables définies dans des boucles
ch_rep, element, nb_Uc, ddlFc, ddlUc, Fc, Uc = [0] * 7

# Unites :
F = input("\nQuelle est l'unité de mesure de force?\t\t")
L = input("Quelle est l'unité de mesure de longueur?\t")
P = input("Quelle est l'unité de mesure de contrainte?\t")
M = f"{F}*{L}"


# ---------
# Fonctions
# ---------

# Calcul du vecteur des charges equivalentes dues a une charge repartie constante appliquee sur un element Poutre1D
def calculer_feq_poutre1d(q, l_poutre):
    feq = np.array([[q * l_poutre / 2],
                    [q * l_poutre ** 2 / 12],
                    [q * l_poutre / 2],
                    [-q * l_poutre ** 2 / 12]])
    return feq


# Calcul par interpolation du deplacement transversal dans un element Poutre1D
def calculer_deplacement_poutre1d(u_tot, ddl, l_poutre, x):
    vi = u_tot[ddl[0] - 1][0]
    ti = u_tot[ddl[1] - 1][0]
    vj = u_tot[ddl[2] - 1][0]
    tj = u_tot[ddl[3] - 1][0]
    dep = ((2 * x ** 3 / l_poutre ** 3 - 3 * x ** 2 / l_poutre ** 2 + 1) * vi +
           (x ** 3 / l_poutre ** 2 - 2 * x ** 2 / l_poutre + x) * ti +
           (-2 * x ** 3 / l_poutre ** 3 + 3 * x ** 2 / l_poutre ** 2) * vj +
           (x ** 3 / l_poutre ** 2 - x ** 2 / l_poutre) * tj)
    return dep


# Calcul de la contrainte dans un element Poutre1D
def calculer_contrainte_poutre1d(u_tot, ddl, e, l_poutre, x, y):
    vi = u_tot[ddl[0] - 1][0]
    ti = Utot[ddl[1] - 1][0]
    vj = Utot[ddl[2] - 1][0]
    tj = Utot[ddl[3] - 1][0]
    sigma = -y * e * ((12 * x / l_poutre ** 3 - 6 / l_poutre ** 2) * vi +
                      (6 * x / l_poutre ** 2 - 4 / l_poutre) * ti +
                      (-12 * x / l_poutre ** 3 + 6 / l_poutre ** 2) * vj +
                      (6 * x / l_poutre ** 2 - 2 / l_poutre) * tj)
    return sigma


# Calcul de l'inertie en Iz pour les calculs
def calcul_Iz():
    colors = {
        'rectangle': '\033[91m',  # Red
        'triangle': '\033[92m',  # Green
        'cercle': '\033[93m',  # Yellow
        'demi-cercle': '\033[94m',  # Blue
        'cercle-mince': '\033[95m',
        'reset': '\033[0m'}  # Reset to default color}  # Magenta
    while True:
        try:
            print(f"{colors['rectangle']}rectangle{colors['reset']}, "
                  f"{colors['triangle']}triangle{colors['reset']}, "
                  f"{colors['cercle']}cercle{colors['reset']}, "
                  f"{colors['demi-cercle']}demi-cercle{colors['reset']}, "
                  f"{colors['cercle-mince']}cercle-mince{colors['reset']}"
                  "\tNone")

            type_element = input("Quelle type de poutre avez-vous :\t")

        except (ValueError, SyntaxError, TypeError):
            print('Erreur de saisie, veuillez recommencer.')
            continue

        if type_element == "rectangle":
            b = float(input(f"Quelle est la valeur de la base en {L} ? "))
            h = float(input(f"Quelle est la valeur de la hauteur en {L} ? "))
            Iz = (b * (h ** 3)) / 12
            print("Iz =", Iz)
            return Iz

        elif type_element == "triangle":
            b = float(input(f"Quelle est la valeur de la base en {L} ? "))
            h = float(input(f"Quelle est la valeur de la hauteur en {L} ? "))
            Iz = (b * (h ** 3)) / 36
            print("Iz =", Iz)
            return Iz

        elif type_element == "cercle":
            r = float(input(f"Quelle est la valeur du rayon en {L} ? "))
            Iz = (math.pi * (r ** 4)) / 4
            print("Iz =", Iz)
            return Iz

        elif type_element == "demi-cercle":
            r = float(input(f"Quelle est la valeur du rayon en {L} ? "))
            Iz = math.pi * (r ** 4) * ((1 / 8) - (8 / 9 * math.pi ** 2))
            print("Iz =", Iz)
            return Iz

        elif type_element == "cercle-mince":
            rm = float(input(f"Quelle est la valeur du rayon moyen en {L} ? "))
            t = float(input(f"Quelle est la valeur de l'épaisseur en {L} ? "))
            Iz = math.pi * (rm ** 3) * t
            print("Iz =", Iz)
            return Iz
        elif type_element == "None":
            Iz = float(input("Valeur de Yz:"))
            print("Iz =", Iz)
            return Iz
        else:
            print('Erreur de saisie, veuillez recommencer.')
            continue


# ----------------------------
# Proprietes de chaque element
# ----------------------------

nb_element = int(input("\nCombien d'éléments contient la structure?\t"))

# Creation d'un dictionnaire avez toutes les cases pour chaque noeuds
elements = {'ddl': [[0, 0, 0, 0] for _ in range(nb_element)], 'L': [0] * nb_element, 'Iz': [0] * nb_element,
            'E': [0] * nb_element, 'q': [0] * nb_element, 'ymax': [0] * nb_element, 'k': [0] * nb_element,
            'feq': [0] * nb_element}

redo = True
while redo:
    for element in range(nb_element):
        print("\n")
        # Request the degrees of freedom of the element
        while True:
            try:
                for i in range(len(elements['ddl'][element])):
                    elements['ddl'][element][i] = int(
                    input(f"Quel est le degré de liberté {i + 1} pour l'élément {element + 1}?\t"))
                break  # Exit this loop after successful input
            except (ValueError, SyntaxError, TypeError):
                print("Entrée invalide, veuillez réessayer.")

        # Request the length of the element
        while True:
            try:
                elements['L'][element] = float(input(f"Quelle est la longueur de l'élément {element + 1} en {L}?\t"))
                break
            except (ValueError, SyntaxError, TypeError):
                print("Entrée invalide, veuillez réessayer.")

        # Calculate Iz
        elements['Iz'][element] = calcul_Iz()

        # Request the modulus of elasticity
        while True:
            try:
                elements['E'][element] = float(
                    input(f"Quel est le module d'élasticité en {P} de l'élément {element + 1}? "))
                break
            except (ValueError, SyntaxError, TypeError):
                print("Entrée invalide, veuillez réessayer.")
        while True:
            try:
                input_ch_rep = input(f"Entrez '1' si l'élément {element + 1}"
                                     f" a une charge répartie, '0' s'il n'en a pas.\t")
                if input_ch_rep == '1' or input_ch_rep == '0':
                    ch_rep = bool(int(input_ch_rep))
                    break
                else:
                    print("Entrée invalide, veuillez réessayer")
                    continue
            except (ValueError, SyntaxError, TypeError):
                print("Entrée invalide, veuillez réessayer")
                continue
        while True:
            if ch_rep:
                q_present = 1
                while True:
                    try:
                        elements['q'][element] = float(
                            input(f"Quelle est la charge répartie en {F}/{L} de l'élément {element + 1}? "))
                        break
                    except (ValueError, SyntaxError, TypeError):
                        print("Entrée invalide, veuillez réessayer.")
                break
            else:
                q_present = 0
                print("Cet élément n'a pas de charge répartie")
                elements['q'][element] = False
                break

        # Request the maximum distance from the neutral axis
        while True:
            try:
                elements['ymax'][element] = float(input(f"Quel est le ymax en {L} de l'élément {element + 1}? "))
                break
            except (ValueError, SyntaxError, TypeError):
                print("Entrée invalide, veuillez réessayer.")

        # Calculate k for the element
        elements['k'][element] = calculer_k_poutre1d(elements['E'][element], elements['Iz'][element],
                                                     elements['L'][element])

        # If distributed load is present, calculate feq
        if q_present == 1:
            elements['feq'][element] = calculer_feq_poutre1d(elements['q'][element], elements['L'][element])
        else:
            elements['feq'][element] = False

    # Check if the user wants to redo the process
    redo_input = input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer: ')
    redo = redo_input.strip() == '1'

# ----------
# Assemblage
# ----------
# Le plus 1 donne le nombre de noeud et *2 car 2 degres de liberte par noeud
Ktot = np.zeros(((nb_element + 1) * 2, ((nb_element + 1) * 2)))
for i in range(nb_element):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])

Feqtot = np.zeros(((nb_element + 1) * 2, 1))
for i in range(nb_element):
    if elements['feq'][i] is False:
        continue
    else:
        Feqtot = assembler_vecteur(Feqtot, elements['feq'][i], elements['ddl'][i])
# -------------------------
# Conditions aux frontieres
# -------------------------
redo = True
while redo is True:
    while True:
        try:
            nb_Uc = int(input('\nCombien de déplacements sont connus?\t'))
        except (ValueError, SyntaxError, TypeError):
            continue
        break
    ddlUc = [0] * nb_Uc
    Uc = np.zeros((nb_Uc, 1))
    print("\n")
    for i in range(nb_Uc):
        while True:
            try:
                ddlUc[i] = int(input(f'Numéro du ddl connu #{i + 1}: '))
                if (ddlUc[i] - 1) % 2 == 0:
                    unite = L
                else:
                    unite = 'rad'
                Uc[i][0] = eval(input(f'Déplacement en {unite} de U{ddlUc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    print("\n")
    for i in range(nb_Uc):
        if (ddlUc[i] - 1) % 2 == 0:
            unite = L
        else:
            unite = 'rad'
        print(f'U{ddlUc[i]}:\t{Uc[i][0]} {unite}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer'))

redo = True
while redo is True:
    try:
        nb_Fc = int(input('\nCombien de forces/moments sont connus?\t'))
    except (ValueError, SyntaxError, TypeError):
        continue
    ddlFc = [0] * nb_Fc
    Fc = np.zeros((nb_Fc, 1))
    print("\n")
    for i in range(nb_Fc):
        while True:
            try:
                ddlFc[i] = int(input(f'Numéro de la force ou du moment connu #{i + 1}: '))
                if (ddlFc[i] - 1) % 2 == 0:
                    unite = F
                else:
                    unite = M
                Fc[i][0] = eval(input(f'Grandeur en {unite} de  F{ddlFc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    print("\n")
    for i in range(nb_Fc):
        if (ddlFc[i] - 1) % 2 == 0:
            unite = F
        else:
            unite = M
        print(f'F{ddlFc[i]}:\t{Fc[i][0]:.2} {unite}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer'))

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Ktot, ddlFc, ddlFc)
Kcc = extraire_matrice(Ktot, ddlFc, ddlUc)
Kii = extraire_matrice(Ktot, ddlUc, ddlFc)
Kci = extraire_matrice(Ktot, ddlUc, ddlUc)
Feqi = extraire_vecteur(Feqtot, ddlUc)
Feqc = extraire_vecteur(Feqtot, ddlFc)

print('\nKtot:')
for i in range(len(Ktot)):
    print(Ktot[i])
print('\nKic:')
for i in range(len(Kic)):
    print(Kic[i])
print('\nKcc:')
for i in range(len(Kcc)):
    print(Kcc[i])
print('\nKii:')
for i in range(len(Kii)):
    print(Kii[i])
print('\nKci:')
for i in range(len(Kci)):
    print(Kci[i])

# --------
# Solution
# --------

Ui = np.linalg.inv(Kic) @ (Fc + Feqc - Kcc @ Uc)
Fi = Kii @ Ui + Kci @ Uc - Feqi

# --------------
# Reconstruction
# --------------

Utot = reconstruire_vecteur(Uc, ddlUc, Ui, ddlFc)
Ftot = reconstruire_vecteur(Fc, ddlFc, Fi, ddlUc)

# --------
# Reponses
# --------

# Print de tous les deplacement
print("\nUtot\n")
for i in range(len(Utot)):
    if i % 2 == 0:
        print(f"U{i + 1} : {Utot[i][0]} {L}")
    else:
        print(f"U{i + 1} : {Utot[i][0]} rad")

# Print de toutes les forces
print("\nFtot\n")
for i in range(len(Ftot)):
    if i % 2 == 0:
        print(f"F{i + 1} : {Ftot[i][0]} {F}")
    else:
        print(f"F{i + 1} : {Ftot[i][0]} {M}")

# Demande a l'utilisateur si il veux des deplacements particulier

while True:
    while True:
        try:
            user_input = input(
                "\nCherchez-vous un déplacement en particulier? Tapez 1 pour oui et 0 pour non :\t")
            if user_input == "1":
                break
            elif user_input == "0":
                redo = False
                break
            else:
                print("Entrée invalide, veuillez entrer 1 ou 0.")
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    if user_input == "0":
        break
    while True:
        try:
            Nddl = int(input("Dans quel élément se trouve le déplacement inconnu?\t"))
            break
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    while True:
        try:
            x = float(input(f"À quelle distance du début de l'élément {Nddl} se trouve le déplacement inconnu?\t"))
            break
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    v = calculer_deplacement_poutre1d(Utot, elements['ddl'][Nddl - 1], elements['L'][Nddl - 1], x)
    print(f"V à {x}{L} du début de l'élément {Nddl} = %.5f mm" % v)

# Demande a l'utilisateur si il veux des contraintes particulier
while True:
    while True:
        try:
            user_input = input(
                "\nCherchez-vous une contrainte en particulier? Tapez 1 pour oui et 0 pour non : ")
            if user_input == "1":
                break
            elif user_input == "0":
                break
            else:
                print("Entrée invalide, veuillez entrer 1 ou 0.")
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    if user_input == "0":
        break
    while True:
        try:
            Nddl = int(input("Dans quel élément se trouve la contrainte inconnue?\t"))
            break
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    while True:
        try:
            x = float(input(f"À quelle distance du début de l'élément {Nddl} se trouve la contrainte inconnue?\t"))
            break
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    s = calculer_contrainte_poutre1d(Utot, elements['ddl'][Nddl - 1], elements['E'][Nddl - 1], elements['L'][Nddl - 1],
                                     x, elements['ymax'][Nddl - 1])
    print(f"SigmaMax à {x}{L} du début de l'élément {Nddl} = %.5f {P}" % s)