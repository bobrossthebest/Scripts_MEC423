import numpy as np

from Modules.Fonctions_partagées import (assembler_matrice, extraire_matrice, extraire_vecteur,
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

# Calcul de la matrice ksi en etat plan de contrainte
def calculer_ksi_epc(e, nu):
    ksi = e / (1 - nu ** 2) * np.array([[1, nu, 0],
                                        [nu, 1, 0],
                                        [0, 0, (1 - nu) / 2]])
    return ksi


# Calcul de l'aire d'un element T3
def calculer_aire_t3(xi, yi, xj, yj, xk, yk):
    temp = np.array([[1, xi, yi],
                     [1, xj, yj],
                     [1, xk, yk]])
    aire = 0.5 * np.linalg.det(temp)
    return aire


# Calcul de la matrice de B d'un element T3
def calculer_b_t3(xi, yi, xj, yj, xk, yk):
    aire = calculer_aire_t3(xi, yi, xj, yj, xk, yk)
    b = 1 / (2 * aire) * np.array([[yj - yk, 0, yk - yi, 0, yi - yj, 0],
                                   [0, xk - xj, 0, xi - xk, 0, xj - xi],
                                   [xk - xj, yj - yk, xi - xk, yk - yi, xj - xi, yi - yj]])
    return b


# Calcul de la matrice de rigidite d'un element T3
def calculer_k_t3(e, nu, t, xi, yi, xj, yj, xk, yk):
    ksi = calculer_ksi_epc(e, nu)
    a = calculer_aire_t3(xi, yi, xj, yj, xk, yk)
    b = calculer_b_t3(xi, yi, xj, yj, xk, yk)
    k = t * a * (np.transpose(b) @ ksi @ b)
    return k


# Evaluation de la fonction d'interpolation Ni de l'element T3 a un point donne
def evaluer_ni_t3(xi, yi, xj, yj, xk, yk, x, y):
    aire = calculer_aire_t3(xi, yi, xj, yj, xk, yk)
    ni = 1 / (2 * aire) * ((xj * yk - xk * yj) + (yj - yk) * x + (xk - xj) * y)
    return ni


# Evaluation de la fonction d'interpolation Nj de l'element T3 a un point donne
def evaluer_nj_t3(xi, yi, xj, yj, xk, yk, x, y):
    aire = calculer_aire_t3(xi, yi, xj, yj, xk, yk)
    nj = 1 / (2 * aire) * ((xk * yi - xi * yk) + (yk - yi) * x + (xi - xk) * y)
    return nj


# Evaluation de la fonction d'interpolation Nk de l'element T3 a un point donne
def evaluer_nk_t3(xi, yi, xj, yj, xk, yk, x, y):
    aire = calculer_aire_t3(xi, yi, xj, yj, xk, yk)
    nk = 1 / (2 * aire) * ((xi * yj - xj * yi) + (yi - yj) * x + (xj - xi) * y)
    return nk


# Calcul des contraintes dans un element T3
def calculer_contraintes_t3(utot, ddl, e, nu, xi, yi, xj, yj, xk, yk):
    ksi = calculer_ksi_epc(e, nu)
    b = calculer_b_t3(xi, yi, xj, yj, xk, yk)
    u = extraire_vecteur(utot, ddl)
    sig = ksi @ b @ u
    return sig


# Calcul de la contrainte equivalente de Mises en etat plan de contrainte
def calculer_mises_epc(sig):
    sx = sig[0][0]
    sy = sig[1][0]
    txy = sig[2][0]
    mises = (sx ** 2 + sy ** 2 + 3 * txy ** 2 - sx * sy) ** 0.5
    return mises


# ----------------------------
# Boucle pour l'entrée des propriétés des noeuds
# ----------------------------


redo = True
while redo is True:
    try:
        nb_noeuds = int(input("\nCombien de noeuds contient la surface?\t"))
    except (ValueError, SyntaxError, TypeError):
        continue
    noeuds = {'x': [0] * nb_noeuds, 'y': [0] * nb_noeuds,
              'ddlx': [0] * nb_noeuds, 'ddly': [0] * nb_noeuds}
    for i in range(nb_noeuds):
        # boucle for pour passer à travers les noeuds, boucle while pour valider les entrées
        while True:
            try:
                noeuds['x'][i] = float(input(f"\nPosition x du noeud {i + 1} en {L}:\t"))
                # commence à 1 quand i est à 0
                noeuds['ddlx'][i] = 2 * i + 1

                noeuds['y'][i] = float(input(f"Position y du noeud {i + 1} en {L}:\t"))
                # commence à 2 quand i est à 0
                noeuds['ddly'][i] = 2 * i + 2

            except (ValueError, SyntaxError, TypeError):
                print('Entrée invalide pour ce noeud, corrigez les valeurs\n')
                continue

            # Si aucun problème d'entrée, sortir de la boucle d'entrée
            break

    print('\n')
    for i in range(nb_noeuds):
        # opérateur ">" pour aligner à droite
        print(f"Noeud {i + 1}:\t[{noeuds['x'][i]:>6}],\t[{noeuds['y'][i]:>6}], "
              f"ddl {noeuds['ddlx'][i]} et {noeuds['ddly'][i]}")
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

# ----------------------------
# Boucle pour l'entrée des propriétés des éléments
# ----------------------------


redo = True
while redo is True:
    try:
        nb_elements = int(input("Combien d'éléments T3 contient la surface? "))
    except (ValueError, SyntaxError, TypeError):
        continue
    vide = [0] * int(nb_elements)
    elements = {'noeud_i': vide.copy(), 'noeud_j': vide.copy(), 'noeud_k': vide.copy(),
                'E': vide.copy(), 't': vide.copy(), 'nu': vide.copy(),
                'ddl': vide.copy(), 'xi': vide.copy(), 'yi': vide.copy(),
                'xj': vide.copy(), 'yj': vide.copy(), 'xk': vide.copy(), 'yk': vide.copy(),
                'A': vide.copy(), 'B': vide.copy(), 'ksi': vide.copy(), 'k': vide.copy()}

    print(f"\nLister les noeuds dans le sens antihoraire autour de l'élément")
    for i in range(nb_elements):
        # soustraction de 1 pour passer du numéro du noeud à son indice dans le tableau
        while True:
            try:
                noeud_i = int(input(f"\n ÉLÉMENT {i + 1}: \nPremier noeud:\t\t")) - 1
                noeud_j = int(input("Second noeud:\t\t")) - 1
                noeud_k = int(input("Troisième noeud:\t")) - 1
                if not (0 <= noeud_i < nb_noeuds and 0 <= noeud_j < nb_noeuds and 0 <= noeud_k < nb_noeuds):
                    print("\t Un des noeuds n'est pas défini.")
                    continue
                # On rajoute 1 pour l'affichage à la fin
                elements['noeud_i'][i] = noeud_i + 1
                elements['noeud_j'][i] = noeud_j + 1
                elements['noeud_k'][i] = noeud_k + 1
            except (ValueError, SyntaxError, TypeError):
                print("Valeur invalide")
                continue

            elements['ddl'][i] = np.array(
                [noeuds['ddlx'][noeud_i], noeuds['ddly'][noeud_i],
                 noeuds['ddlx'][noeud_j], noeuds['ddly'][noeud_j],
                 noeuds['ddlx'][noeud_k], noeuds['ddly'][noeud_k]])
            elements['xi'][i], elements['yi'][i] = (
                noeuds['x'][noeud_i], noeuds['y'][noeud_i])
            elements['xj'][i], elements['yj'][i] = (
                noeuds['x'][noeud_j], noeuds['y'][noeud_j])
            elements['xk'][i], elements['yk'][i] = (
                noeuds['x'][noeud_k], noeuds['y'][noeud_k])

            if calculer_aire_t3(elements['xi'][i], elements['yi'][i],
                                elements['xj'][i], elements['yj'][i],
                                elements['xk'][i], elements['yk'][i]) < 0:
                print("Les noeuds n'étaient pas dans le sens horaire. Leur ordre a été changé.")

                # Changement des ddl de l'élément
                elements['ddl'][i][2], elements['ddl'][i][3], elements['ddl'][i][4], elements['ddl'][i][5] = (
                    elements['ddl'][i][4], elements['ddl'][i][5], elements['ddl'][i][2], elements['ddl'][i][3])

                # Changement des positions de l'élément
                elements['xj'][i], elements['yj'][i], elements['xk'][i], elements['yk'][i] = (
                    elements['xk'][i], elements['yk'][i], elements['xj'][i], elements['yj'][i])

                # Changement des marqueurs de noeuds
                elements['noeud_j'][i], elements['noeud_k'][i] = elements['noeud_k'][i], elements['noeud_j'][i]

                print(f"Le nouvel ordre des noeuds pour cet élément est donc {elements['noeud_i'][i]}, "
                      f"{elements['noeud_j'][i]} et {elements['noeud_k'][i]}.")
            try:
                elements['t'][i] = eval(input(f"Épaisseur en {L}:\t"))
                elements['E'][i] = eval(input(f"Module d'élasticité en {P}:\t"))
                elements['nu'][i] = eval(input(f"Module de Poisson:\t"))
            except (SyntaxError, ValueError, TypeError):
                print("Erreur dans les valeurs entrées")
                continue
            elements['A'][i] = calculer_aire_t3(elements['xi'][i], elements['yi'][i],
                                                elements['xj'][i], elements['yj'][i],
                                                elements['xk'][i], elements['yk'][i])

            elements['B'][i] = calculer_b_t3(elements['xi'][i], elements['yi'][i],
                                             elements['xj'][i], elements['yj'][i],
                                             elements['xk'][i], elements['yk'][i])

            elements['ksi'][i] = calculer_ksi_epc(elements['E'][i], elements['nu'][i])

            elements['k'][i] = calculer_k_t3(elements['E'][i], elements['nu'][i], elements['t'][i],
                                             elements['xi'][i], elements['yi'][i],
                                             elements['xj'][i], elements['yj'][i],
                                             elements['xk'][i], elements['yk'][i])
            break

    print('\n')
    print('#\ti-j-k\t       DDL        \t xi \t yi \t xj \t yj \t xk \t yk \t t \t  E  \t v ')
    for i in range(nb_elements):
        print(f"{i + 1}\t{elements['noeud_i'][i]}-{elements['noeud_j'][i]}-{elements['noeud_k'][i]}\t"
              f"{elements['ddl'][i][0]:<2} {elements['ddl'][i][1]:<2} {elements['ddl'][i][2]:<2} "
              f"{elements['ddl'][i][3]:<2} {elements['ddl'][i][4]:<2} {elements['ddl'][i][5]:<2}\t"
              f"{elements['xi'][i]:<4}\t{elements['yi'][i]:<4}\t{elements['xj'][i]:<4}\t{elements['yj'][i]:<4}\t"
              f"{elements['xk'][i]:<4}\t{elements['yk'][i]:<4}\t{elements['t'][i]:<3}\t{int(elements['E'][i]):<6}\t"
              f"{elements['nu'][i]:<3}")
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

    for i in range(nb_elements):
        print(f'\n\tÉlément {i + 1}:')
        print(f'A: {elements['A'][i]:.0f}')
        print(f'\nB: 1/(2*{elements['A'][i]:.0})*')
        for j in range(len(elements['B'][i])):
            print(f'{2 * elements['A'][i] * elements['B'][i][j]}')
        print(f'\nksi:')
        for k in range(len(elements['ksi'][i][0])):
            print(f'{elements['ksi'][i][k]}')

# ----------
# Assemblage
# ----------

Ktot = np.zeros((nb_noeuds * 2, nb_noeuds * 2))
for i in range(nb_elements):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])

# -------------------------
# Conditions aux frontieres
# -------------------------

redo = True
while redo is True:
    try:
        nb_Uc = int(input('\nCombien de déplacements sont connus? '))
    except (ValueError, SyntaxError, TypeError):
        continue
    ddlUc = [0] * nb_Uc
    Uc = np.zeros((nb_Uc, 1))
    for i in range(nb_Uc):
        while True:
            try:
                ddlUc[i] = int(input(f'Numéro du ddl connu #{i + 1}: '))
                Uc[i][0] = eval(input(f'Déplacement en {L} de U{ddlUc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    print('\n')
    for i in range(nb_Uc):
        print(f'U{ddlUc[i]}:\t{Uc[i][0]} {L}')
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
                Fc[i][0] = eval(input(f'Grandeur en {F} de  F{ddlFc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    print('\n')
    for i in range(nb_Fc):
        print(f'F{ddlFc[i]}:\t{Fc[i][0]:.2} {F}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Ktot, ddlFc, ddlFc)
Kcc = extraire_matrice(Ktot, ddlFc, ddlUc)
Kii = extraire_matrice(Ktot, ddlUc, ddlFc)
Kci = extraire_matrice(Ktot, ddlUc, ddlUc)

print('Ktot:')
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
print('\n')

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

# Déplacements interpolés
while True:
    while True:
        try:
            user_input = input(
                "\nCherchez-vous un déplacement spécifique? Tapez 1 pour oui et 0 pour non :\t")
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
            # -1 pour passer du numérotage à l'indice dans le tableau
            N = int(input("Dans quel élément se trouve le déplacement recherché?\t")) - 1
            break
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    while True:
        try:
            x = float(input(f"À quelle position x (en {L}) se trouve le déplacement inconnu?\t"))
            y = float(input(f"À quelle position y (en {L})?\t"))
            break
        except (ValueError, SyntaxError, TypeError):
            print("Entrée invalide, veuillez réessayer.")
    Ni = evaluer_ni_t3(elements['xi'][N], elements['yi'][N], elements['xj'][N], elements['yj'][N],
                       elements['xk'][N], elements['yk'][N], x, y)
    Nj = evaluer_nj_t3(elements['xi'][N], elements['yi'][N], elements['xj'][N], elements['yj'][N],
                       elements['xk'][N], elements['yk'][N], x, y)
    Nk = evaluer_nk_t3(elements['xi'][N], elements['yi'][N], elements['xj'][N], elements['yj'][N],
                       elements['xk'][N], elements['yk'][N], x, y)

    u = [Utot[elements['ddl'][N][0] - 1][0], Utot[elements['ddl'][N][2] - 1][0], Utot[elements['ddl'][N][4] - 1][0]]
    v = [Utot[elements['ddl'][N][1] - 1][0], Utot[elements['ddl'][N][3] - 1][0], Utot[elements['ddl'][N][5] - 1][0]]

    Nijk = [Ni, Nj, Nk]
    up = np.dot(u, Nijk)
    vp = np.dot(v, Nijk)

    print(f"Déplacements au point {x}, {y} dans l'élément {N+1} : U = {up} {L}, V = {vp} {L}")

# Tableau des contraintes, incluant Von Mises
tab_sigma = []
for i in range(nb_elements):
    tab_sigma.append([0, 0])
for i in range(nb_elements):
    tab_sigma[i][0] = calculer_contraintes_t3(Utot, elements['ddl'][i], elements['E'][i], elements['nu'][i],
                                              elements['xi'][i], elements['yi'][i],
                                              elements['xj'][i], elements['yj'][i],
                                              elements['xk'][i], elements['yk'][i])
    tab_sigma[i][1] = calculer_mises_epc(tab_sigma[i][0])

print('\nContraintes:\n')
for i in range(nb_elements):
    print(f"Élément {i+1}: \u03c3x = {tab_sigma[i][0][0][0]:.3},\t\u03c3y = {tab_sigma[i][0][1][0]:.3},\t"
          f"\U0001d70fxy = {tab_sigma[i][0][2][0]:.3}")
    print(f"\t\t\u03c3 de von Mises: {tab_sigma[i][1]:.3}")
    print(f"\t\tFS = {25 / tab_sigma[i][1]:.3}\n")
