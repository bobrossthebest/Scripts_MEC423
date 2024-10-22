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

    print(f"\tLister les noeuds dans le sens antihoraire autour de l'élément")
    for i in range(nb_elements):
        # soustraction de 1 pour passer du numéro du noeud à son indice dans le tableau
        while True:
            try:
                noeud_i = int(input(f"\n ÉLÉMENT {i + 1}: \nPremier noeud:\t")) - 1
                noeud_j = int(input("Second noeud:\t")) - 1
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
                elements['nu'][i] = eval(input(f"Module de Poisson:\t\t"))
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
    print('#\ti-j-k\t       DDL        \txi\tyi\txj\tyj\txk\tyk\t t \t E  \t v ')
    for i in range(nb_elements):
        print(f"{i + 1}\t{elements['noeud_i'][i]}-{elements['noeud_j'][i]}-{elements['noeud_k'][i]}\t"
              f"{float(elements['ddl'][i][0]):<2} {float(elements['ddl'][i][1]):<2} {float(elements['ddl'][i][2]):<2} "
              f"{float(elements['ddl'][i][3]):<2} {float(elements['ddl'][i][4]):<2} {float(elements['ddl'][i][5]):<2} "
              f"{elements['xi'][i]:<2}\t{elements['yi'][i]:<2}\t{elements['xj'][i]:<2}\t{elements['xj'][i]:<2}\t"
              f"{elements['xk'][i]:<2}\t{elements['yk'][i]:<2}\t{float(elements['t'][i]):<3f}\t {elements['nu'][i]:<3}")
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1 = np.array([1, 2, 5, 6, 3, 4])
xi1, yi1 = 0, 0
xj1, yj1 = 20, 0
xk1, yk1 = 0, 10
t1 = 1.2
E1 = 2e5
nu1 = 0.3
k1 = calculer_k_t3(E1, nu1, t1, xi1, yi1, xj1, yj1, xk1, yk1)

ddl2 = np.array([3, 4, 5, 6, 7, 8])
xi2, yi2 = 0, 10
xj2, yj2 = 20, 0
xk2, yk2 = 20, 10
t2 = 1.2
E2 = 2e5
nu2 = 0.3
k2 = calculer_k_t3(E2, nu2, t2, xi2, yi2, xj2, yj2, xk2, yk2)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((8, 8))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([1, 2, 3, 4])
Uc = np.array([[0], [0], [0], [0]])
ddlFc = np.array([5, 6, 7, 8])
Fc = np.array([[0], [0], [0], [-3000]])

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Ktot, ddlFc, ddlFc)
Kcc = extraire_matrice(Ktot, ddlFc, ddlUc)
Kii = extraire_matrice(Ktot, ddlUc, ddlFc)
Kci = extraire_matrice(Ktot, ddlUc, ddlUc)

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

xP = 17
yP = 8
Ni2P = evaluer_ni_t3(xi2, yi2, xj2, yj2, xk2, yk2, xP, yP)
Nj2P = evaluer_nj_t3(xi2, yi2, xj2, yj2, xk2, yk2, xP, yP)
Nk2P = evaluer_nk_t3(xi2, yi2, xj2, yj2, xk2, yk2, xP, yP)
vi2 = extraire_vecteur(Utot, [4])
vj2 = extraire_vecteur(Utot, [6])
vk2 = extraire_vecteur(Utot, [8])
vP = vi2[0][0] * Ni2P + vj2[0][0] * Nj2P + vk2[0][0] * Nk2P
print('vP = %.4f mm' % vP)

sig1 = calculer_contraintes_t3(Utot, ddl1, E1, nu1, xi1, yi1, xj1, yj1, xk1, yk1)
mises1 = calculer_mises_epc(sig1)
print('Mises_1 = %.1f MPa' % mises1)

sig2 = calculer_contraintes_t3(Utot, ddl2, E2, nu2, xi2, yi2, xj2, yj2, xk2, yk2)
mises2 = calculer_mises_epc(sig2)
print('Mises_2 = %.1f MPa' % mises2)
