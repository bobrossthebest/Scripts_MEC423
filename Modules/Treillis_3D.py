import numpy as np

from Modules.Fonctions_partagées import (assembler_matrice, extraire_matrice,
                                         reconstruire_vecteur)

# Unites :
F = input("\nQuelle est l'unité de mesure de force?\t\t")
L = input("Quelle est l'unité de mesure de longueur?\t")
P = input("Quelle est l'unité de mesure de contrainte?\t")

# Variables définies dans des boucles
noeuds, nb_noeuds, nb_elements, elements, ddlFc, ddlUc, Fc, Uc = [False] * 8


# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Barre3D
def calculer_k_barre3d(e, a, xi, yi, zi, xj, yj, zj):
    l_ = ((xj - xi) ** 2 + (yj - yi) ** 2 + (zj - zi) ** 2) ** 0.5
    cx = (xj - xi) / l_
    cy = (yj - yi) / l_
    cz = (zj - zi) / l_
    k = (e * a / l_) * np.array([[cx ** 2, cx * cy, cx * cz, -cx ** 2, -cx * cy, -cx * cz],
                                 [cx * cy, cy ** 2, cy * cz, -cx * cy, -cy ** 2, -cy * cz],
                                 [cx * cz, cy * cz, cz ** 2, -cx * cz, -cy * cz, -cz ** 2],
                                 [-cx ** 2, -cx * cy, -cx * cz, cx ** 2, cx * cy, cx * cz],
                                 [-cx * cy, -cy ** 2, -cy * cz, cx * cy, cy ** 2, cy * cz],
                                 [-cx * cz, -cy * cz, -cz ** 2, cx * cz, cy * cz, cz ** 2]])
    return k


# Calcul de la contrainte dans un element Barre3D
def calculer_contrainte_barre3d(u_tot, ddl, e, alpha, dt, xi, yi, zi, xj, yj, zj):
    ui = u_tot[ddl[0] - 1][0]
    vi = u_tot[ddl[1] - 1][0]
    wi = u_tot[ddl[2] - 1][0]
    uj = u_tot[ddl[3] - 1][0]
    vj = u_tot[ddl[4] - 1][0]
    wj = u_tot[ddl[5] - 1][0]
    l_ = ((xj - xi) ** 2 + (yj - yi) ** 2 + (zj - zi) ** 2) ** 0.5
    cx = (xj - xi) / l_
    cy = (yj - yi) / l_
    cz = (zj - zi) / l_
    sigma = (e / l_) * np.array([[-cx, -cy, -cz, cx, cy, cz]]) @ np.array(
        [[ui], [vi], [wi], [uj], [vj], [wj]]) - np.array([[e * alpha * dt]])
    return sigma[0][0]


# ----------------------------
# Proprietes de chaque element
# ----------------------------

# Boucle d'entrée pour les noeuds

redo = True
while redo is True:
    try:
        nb_noeuds = int(input("\nCombien de noeuds contient la structure?\t"))
    except (ValueError, SyntaxError, TypeError):
        continue
    noeuds = {'x': [0] * nb_noeuds, 'y': [0] * nb_noeuds, 'z': [0] * nb_noeuds,
              'ddlx': [0] * nb_noeuds, 'ddly': [0] * nb_noeuds, 'ddlz': [0] * nb_noeuds, }
    for i in range(nb_noeuds):
        # boucle for pour passer à travers les noeuds, boucle while pour valider les entrées
        while True:
            try:
                noeuds['x'][i] = float(input(f"\nPosition x du noeud {i + 1} en {L}:\t"))
                # commence à 1 quand i est à 0
                noeuds['ddlx'][i] = 3 * i + 1

                noeuds['y'][i] = float(input(f"Position y du noeud {i + 1} en {L}:\t"))
                # commence à 2 quand i est à 0
                noeuds['ddly'][i] = 3 * i + 2

                noeuds['z'][i] = float(input(f"Position z du noeud {i + 1} en {L}:\t"))
                # commence à 3 quand i est à 0
                noeuds['ddlz'][i] = 3 * i + 3

            except (ValueError, SyntaxError, TypeError):
                print('Entrée invalide pour ce noeud, corrigez\n')
                continue

            # Si aucun problème d'entrée, sortir de la boucle d'entrée
            break

    print('\n')
    for i in range(nb_noeuds):
        # opérateur ">" pour aligner à droite
        print(f"Noeud {i+1}:\t[{noeuds['x'][i]:>6}],\t[{noeuds['y'][i]:>6}],\t[{noeuds['z'][i]:>6}], "
              f"ddl {noeuds['ddlx'][i]} à {noeuds['ddlz'][i]}")
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

# Boucle d'entrée pour les éléments

redo = True
while redo is True:
    try:
        nb_elements = int(input("Combien d'éléments contient la structure? "))
    except (ValueError, SyntaxError, TypeError):
        continue
    vide = [0] * int(nb_elements)
    elements = {'ddl': vide.copy(), 'xi': vide.copy(), 'yi': vide.copy(), 'zi': vide.copy(),
                'xj': vide.copy(), 'yj': vide.copy(), 'zj': vide.copy(),
                'E': vide.copy(), 'A': vide.copy(), 'alpha': vide.copy(), 'dT': vide.copy(),
                'k': vide.copy(), 'feq': vide.copy(), 'noeud_i': vide.copy(), 'noeud_j': vide.copy()}

    print(f"\tPour un ressort, poser un module d'élasticité de 0")
    for i in range(nb_elements):
        # soustraction de 1 pour passer du numéro du noeud à son indice dans le tableau
        while True:
            try:
                noeud_i = int(input(f"\n ÉLÉMENT {i + 1}: \nNoeud avant l'élément:\t")) - 1
                noeud_j = int(input("Noeud après l'élément:\t")) - 1
                elements['noeud_i'][i] = noeud_i
                elements['noeud_j'][i] = noeud_j
            except (ValueError, SyntaxError, TypeError):
                print("Valeur invalide")
                continue

            elements['ddl'][i] = np.array(
                [noeuds['ddlx'][noeud_i], noeuds['ddly'][noeud_i], noeuds['ddlz'][noeud_i],
                 noeuds['ddlx'][noeud_j], noeuds['ddly'][noeud_j], noeuds['ddlz'][noeud_j]])
            elements['xi'][i], elements['yi'][i], elements['zi'][i] = (
                noeuds['x'][noeud_i], noeuds['y'][noeud_i], noeuds['z'][noeud_i])
            elements['xj'][i], elements['yj'][i], elements['zj'][i] = (
                noeuds['x'][noeud_j], noeuds['y'][noeud_j], noeuds['z'][noeud_j])

            try:
                elements['E'][i] = eval(input(f"Module d'élasticité en {P}:\t"))
                if elements['E'][i] > 0:
                    elements['A'][i] = eval(input(f"Aire de section en {L}^2:\t"))
                else:
                    elements['A'][i] = float(input(f"Raideur du ressort en {F}/{L}:\t"))
                elements['dT'][i] = float(input('Différence de température:\t'))
                if elements['dT'][i] != 0:
                    elements['alpha'][i] = float(input("Coefficient de dilatation thermique:\t"))
            except (SyntaxError, ValueError, TypeError):
                print("Erreur dans les valeurs entrées")
                continue
            elements['k'][i] = calculer_k_barre3d(elements['E'][i], elements['A'][i],
                                                  elements['xi'][i], elements['yi'][i], elements['zi'][i],
                                                  elements['xj'][i], elements['yj'][i], elements['zj'][i])
            break

    print('\n')
    for i in range(nb_elements):
        if elements['E'][i] > 0 and elements['dT'][i] == 0:
            print(f"Élément {i + 1} du noeud {elements['noeud_i'][i]} au noeud {elements['noeud_j'][i]}:\t"
                  f"E = {elements['E'][i]:<5} {P},\tA = {elements['A'][i]:<5} {L}^2")
        elif elements['dT'][i] != 0:
            print(f"Élément {i + 1} du noeud {elements['noeud_i'][i]} au noeud {elements['noeud_j'][i]}:\t"
                  f"E = {elements['E'][i]:<5} {P},\tA = {elements['A'][i]:<5} {L}^2,\t"
                  f"dT = {elements['dT'][i]:<5},\talpha = {elements['alpha'][i]:<5}")
        elif elements['E'][i] == 0:
            print(f"Élément {i + 1} du noeud {elements['noeud_i'][i]} au noeud {elements['noeud_j'][i]}:\t"
                  f"k = {elements['A'][i]:<5} {F}/{L}")
        else:
            print(f"Élément {i + 1} du noeud {elements['noeud_i'][i]} au noeud {elements['noeud_j'][i]}:\t"
                  f"k = {elements['A'][i]:<5} {F}/{L},\t"
                  f"dT = {elements['dT'][i]:<5},\talpha = {elements['alpha'][i]:<5}")
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

# ----------
# Assemblage
# ----------

# 3 pcq 3 ddl par noeud
Ktot = np.zeros((nb_noeuds * 3, nb_noeuds * 3))
for i in range(nb_elements):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])

# -------------------------
# Conditions aux frontieres
# -------------------------

redo = True
while redo is True:
    try:
        nb_Uc = int(input('Combien de déplacements sont connus? '))
    except (ValueError, SyntaxError, TypeError):
        continue
    ddlUc = [0] * nb_Uc
    Uc = np.zeros((nb_Uc, 1))
    for i in range(nb_Uc):
        while True:
            try:
                ddlUc[i] = int(input(f'Numéro du ddl connu #{i + 1}: '))
                Uc[i][0] = eval(input(f'Déplacement en {L} du noeud {ddlUc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    for i in range(nb_Uc):
        print(f'Noeud {ddlUc[i]}:\t{Uc[i][0]} {L}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

redo = True
while redo is True:
    try:
        nb_Fc = int(input('Combien de forces sont connues? '))
    except (ValueError, SyntaxError, TypeError):
        continue
    ddlFc = [0] * nb_Fc
    Fc = np.zeros((nb_Fc, 1))
    for i in range(nb_Fc):
        while True:
            try:
                ddlFc[i] = int(input(f'Numéro de la force connue #{i + 1}: '))
                Fc[i][0] = eval(input(f'Grandeur en {F} de la force {ddlFc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    for i in range(nb_Fc):
        print(f'Noeud {ddlFc[i]}:\t{Fc[i][0]:.2} {F}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Ktot, ddlFc, ddlFc)
Kcc = extraire_matrice(Ktot, ddlFc, ddlUc)
Kii = extraire_matrice(Ktot, ddlUc, ddlFc)
Kci = extraire_matrice(Ktot, ddlUc, ddlUc)

print('\n')
for tableau in (Kic, Kcc, Kii, Kci):
    for i in range(len(tableau)):
        for j in range(len(tableau[0])):
            print(tableau[i][j])

# --------
# Solution
# --------

Ui = np.linalg.inv(Kic) @ (Fc - Kcc @ Uc)
Fi = Kii @ Ui + Kci @ Uc

# --------------
# Reconstruction
# --------------

Utot = reconstruire_vecteur(Uc, ddlUc, Ui, ddlFc)
Ftot = reconstruire_vecteur(Fc, ddlFc, Fi, ddlUc)

# --------
# Reponses
# --------

# Tableau de contraintes
tab_sigma = [0] * nb_elements
for i in range(nb_elements):
    tab_sigma[i] = calculer_contrainte_barre3d(Utot, elements['ddl'][i], elements['E'][i],
                                               elements['alpha'][i], elements['dT'][i],
                                               elements['xi'][i], elements['yi'][i], elements['zi'][i],
                                               elements['xj'][i], elements['yj'][i], elements['zj'][i])
# Tableau de forces
tab_force = [0] * nb_elements
for i in range(nb_elements):
    tab_force[i] = elements['A'][i] * tab_sigma[i]

for i in range(nb_elements):
    print(f"Élément {i + 1}\tForce: {tab_force[i]}{F}\tContrainte: {tab_sigma[i]}{P}")
print('\n')
for i in range(len(Ui)):
    print(f"Déplacement {ddlFc[i]}:\t{Ui[i]}")
