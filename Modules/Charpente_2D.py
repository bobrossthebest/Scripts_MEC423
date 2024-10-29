# Unites : N-mm-MPa

import numpy as np
from Modules.Fonctions_partagées import (assembler_matrice, extraire_matrice, extraire_vecteur, reconstruire_vecteur,
                                         calculer_k_poutre2d, calcul_Iz, calcul_q)

# Unites :
F = input("\nQuelle est l'unité de mesure de force?\t\t")
L = input("Quelle est l'unité de mesure de longueur?\t")
P = input("Quelle est l'unité de mesure de contrainte?\t")
M = f"{F}*{L}"

# ---------
# Fonctions
# ---------

# Calcul du vecteur des charges equivalentes dues a une charge repartie constante appliquee sur un element Poutre2D
def calculer_feq_poutre2d(e, a, alpha, d_t, qx, qy, xi, yi, xj, yj):
    l_ = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / l_
    cy = (yj - yi) / l_
    feq_t = a * e * alpha * d_t * np.array([[-cx],
                                           [-cy],
                                           [0],
                                           [cx],
                                           [cy],
                                           [0]])
    feq_q = np.array([[qx * l_ / 2],
                     [qy * l_ / 2],
                     [(qy * cx - qx * cy) * l_ ** 2 / 12],
                     [qx * l_ / 2],
                     [qy * l_ / 2],
                     [(qx * cy - qy * cx) * l_ ** 2 / 12]])
    feq = feq_t + feq_q
    return feq


# Assemblage d'un vecteur
def assembler_vecteur(ftot, f, lig):
    for i in range(len(lig)):
        ftot[lig[i] - 1][0] += f[i][0]
    return ftot


# Calcul des contraintes aux extremites d'un element Poutre2D
def calculer_contraintes_poutre2d(utot, ddl, e, alpha, d_t, xi, yi, xj, yj, yplus, ymoins):
    l_ = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / l_
    cy = (yj - yi) / l_
    r = np.array([[cx, cy, 0, 0, 0, 0],
                  [-cy, cx, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, cx, cy, 0],
                  [0, 0, 0, -cy, cx, 0],
                  [0, 0, 0, 0, 0, 1]])
    disp_global = np.array([[utot[ddl[0] - 1][0]],
                           [utot[ddl[1] - 1][0]],
                           [utot[ddl[2] - 1][0]],
                           [utot[ddl[3] - 1][0]],
                           [utot[ddl[4] - 1][0]],
                           [utot[ddl[5] - 1][0]]])
    disp_local = r @ disp_global
    ui = disp_local[0][0]
    vi = disp_local[1][0]
    ti = disp_local[2][0]
    uj = disp_local[3][0]
    vj = disp_local[4][0]
    tj = disp_local[5][0]
    sig_trac = e * (uj - ui) / l_ - e * alpha * d_t
    sig_flexi_plus = -yplus * e * (6 * (vj - vi) / l_ ** 2 - 2 * (2 * ti + tj) / l_)
    sig_flexi_moins = -ymoins * e * (6 * (vj - vi) / l_ ** 2 - 2 * (2 * ti + tj) / l_)
    sig_flexj_plus = -yplus * e * (6 * (vi - vj) / l_ ** 2 + 2 * (ti + 2 * tj) / l_)
    sig_flexj_moins = -ymoins * e * (6 * (vi - vj) / l_ ** 2 + 2 * (ti + 2 * tj) / l_)
    sig_combi_plus = sig_flexi_plus + sig_trac
    sig_combi_moins = sig_flexi_moins + sig_trac
    sig_combj_plus = sig_flexj_plus + sig_trac
    sig_combj_moins = sig_flexj_moins + sig_trac
    sig = np.array(
        [sig_trac, sig_flexi_plus, sig_flexi_moins, sig_flexj_plus, sig_flexj_moins,
         sig_combi_plus, sig_combi_moins, sig_combj_plus, sig_combj_moins])
    return sig

# ----------------------------
# Boucle pour l'entrée des propriétés des noeuds
# ----------------------------

redo = True
while redo is True:
    try:
        nb_noeuds = int(input("\nCombien de noeuds contient la charpente?\t"))
    except (ValueError, SyntaxError, TypeError):
        continue
    noeuds = {'x': [0] * nb_noeuds, 'y': [0] * nb_noeuds,
              'ddlx': [0] * nb_noeuds, 'ddly': [0] * nb_noeuds, 'ddltheta':[0]*nb_noeuds}
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

                # commence à 3 quand i est à 0
                noeuds['ddltheta'][i] = 3 * i + 3

            except (ValueError, SyntaxError, TypeError):
                print('Entrée invalide pour ce noeud, corrigez les valeurs\n')
                continue

            # Si aucun problème d'entrée, sortir de la boucle d'entrée
            break

    print('\n')
    for i in range(nb_noeuds):
        # opérateur ">" pour aligner à droite
        print(f"Noeud {i + 1}:\t[{noeuds['x'][i]:>6}],\t[{noeuds['y'][i]:>6}], "
              f"ddl {noeuds['ddlx'][i]}, {noeuds['ddly'][i]}, {noeuds['ddltheta'][i]} ")
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

# ----------------------------
# Proprietes de chaque element
# ----------------------------

redo = True
while redo is True:
    try:
        nb_element = int(input("\nCombien d'éléments contient la structure?\t"))
    except (ValueError, SyntaxError, TypeError):
        continue
    vide = [0] * int(nb_element)
    # Creation d'un dictionnaire avez toutes les cases pour chaque noeuds
    elements = {'ddl': [[0, 0, 0, 0, 0, 0] for _ in range(nb_element)],
                'noeud_i': vide.copy(),
                'noeud_j': vide.copy(),
                'xi': vide.copy(),
                'yi': vide.copy(),
                'xj': vide.copy(),
                'yj': vide.copy(),
                'A' : vide.copy(),
                'Sy': vide.copy(),
                'Iz': vide.copy(),
                'yplus': vide.copy(),
                'ymoins': vide.copy(),
                'E': vide.copy(),
                'qx': vide.copy(),
                'qy': vide.copy(),
                'alpha': vide.copy(),
                'dT': vide.copy(),
                'k': vide.copy(),
                'feq': vide.copy()}

    for i in range(nb_element):
        # soustraction de 1 pour passer du numéro du noeud à son indice dans le tableau
        while True:
            try:
                noeud_i = int(input(f"\n ÉLÉMENT {i + 1}: \nPremier noeud:\t\t")) - 1
                noeud_j = int(input("Second noeud:\t\t")) - 1

                if not (0 <= noeud_i < nb_noeuds and 0 <= noeud_j < nb_noeuds):
                    print("\t Un des noeuds n'est pas défini.")
                    continue
                # On rajoute 1 pour l'affichage à la fin
                elements['noeud_i'][i] = noeud_i + 1
                elements['noeud_j'][i] = noeud_j + 1
                elements['xi'][i] = noeuds['x'][noeud_i]
                elements['yi'][i] = noeuds['y'][noeud_i]
                elements['xj'][i] = noeuds['x'][noeud_j]
                elements['yj'][i] = noeuds['y'][noeud_j]
                elements['ddl'][i] = [noeuds['ddlx'][noeud_i],
                                      noeuds['ddly'][noeud_i],
                                      noeuds['ddltheta'][noeud_i],
                                      noeuds['ddlx'][noeud_j],
                                      noeuds['ddly'][noeud_j],
                                      noeuds['ddltheta'][noeud_j]
                                      ]

            except (ValueError, SyntaxError, TypeError):
                print("Valeur invalide")
                continue
            break
        while True:
            try:
                elements['Sy'][i] = eval(input(f"Limite d'ecoulement en {P}:\t"))
                elements['A'][i] = eval(input(f"Aire de la section en {L}\u00B2 :\t"))
                elements['Iz'][i] = calcul_Iz(L)
            except (SyntaxError, ValueError, TypeError):
                print("Erreur dans les valeurs entrées")
                continue
            break
        while True:
            try:
                elements['yplus'][i] = eval(input(f"Y plus de la section en {L}:\t"))
                elements['ymoins'][i] = eval(input(f"Y moins de la section en {L}:\t"))
                elements['E'][i] = eval(input(f"Module d'élasticité en {P}:\t"))
                elements['qx'][i], elements['qy'][i] = calcul_q(i,F,L)
                elements['dT'][i] = float(input('Différence de température:\t'))
                if elements['dT'][i] != 0:
                    elements['alpha'][i] = eval(input("Coefficient de dilatation thermique:\t"))
                elements['k'][i] = calculer_k_poutre2d(elements['E'][i],
                                                       elements['A'][i],
                                                       elements['Iz'][i],
                                                       elements['xi'][i],
                                                       elements['yi'][i],
                                                       elements['xj'][i],
                                                       elements['yj'][i])
                elements['feq'][i] = calculer_feq_poutre2d(elements['E'][i],
                                                           elements['A'][i],
                                                           elements['alpha'][i],
                                                           elements['dT'][i],
                                                           elements['qx'][i],
                                                           elements['qy'][i],
                                                           elements['xi'][i],
                                                           elements['yi'][i],
                                                           elements['xj'][i],
                                                           elements['yj'][i])
            except (SyntaxError, ValueError, TypeError):
                print("Erreur dans les valeurs entrées")
                continue
            break

    print('\n')
    print('#\ti-j\t       DDL        \t xi \t yi \t xj \t yj \t A \t Iz \t  E  \t Sy \t alpha \t dT \t qx \t qy ')
    for i in range(nb_element):
        print(f"{i + 1}\t{elements['noeud_i'][i]}-{elements['noeud_j'][i]}\t\t"
              f"{elements['ddl'][i][0]:<2} {elements['ddl'][i][1]:<2} {elements['ddl'][i][2]:<2} "
              f"{elements['ddl'][i][3]:<2} {elements['ddl'][i][4]:<2} {elements['ddl'][i][5]:<2}\t"
              f"{elements['xi'][i]:<4}\t{elements['yi'][i]:<4}\t{elements['xj'][i]:<4}\t{elements['yj'][i]:<4}\t"
              f"{elements['A'][i]:<4}\t{elements['Iz'][i]:<4}\t{int(elements['E'][i]):<6}\t{elements['Sy'][i]:<4}\t{elements['alpha'][i]:<4}\t"
              f"{elements['dT'][i]:<4}\t{elements['qx'][i]:<4.3f}\t{elements['qy'][i]:<4.3f}\t")

    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))


# ----------
# Assemblage
# ----------

Ktot = np.zeros((nb_noeuds * 3, nb_noeuds * 3))
for i in range(nb_element):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])

Feqtot = np.zeros((nb_noeuds * 3, 1))
for i in range(nb_element):
    Feqtot = assembler_vecteur(Feqtot,  elements['feq'][i], elements['ddl'][i])


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
        if ddlUc[i] % 3 == 0:  # Check si c'est le c'est un angle ou deplacement
            print(f'U{ddlUc[i]} :\t{Uc[i][0]} \u03B8')
        else:
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
        if ddlFc[i] % 3 == 0:  # Check si c'est le c'est un angle ou deplacement
            print(f'U{ddlFc[i]} :\t{Fc[i][0]} {F}*{L}')
        else:
            print(f'U{ddlFc[i]}:\t{Fc[i][0]} {L}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

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
    if (i+1) % 3 == 0:
        print(f"U{i + 1} : {Utot[i][0]:.3f} \u03B8")
    else:
        print(f"U{i + 1} : {Utot[i][0]:.3f} {L}")

# Print de toutes les forces
print("\nFtot\n")
for i in range(len(Ftot)):
    if (i+1) % 3 == 0:
        print(f"F{i + 1} : {Ftot[i][0]} {M}")
    else:
        print(f"F{i + 1} : {Ftot[i][0]} {F}")

# Tableau des contraintes, incluant Von Mises
tab_sigma_contrainte = []
for i in range(nb_element):
    tab_sigma_contrainte.append([0, 0])
for i in range(nb_element):
    tab_sigma_contrainte[i][0] = max(abs(calculer_contraintes_poutre2d(Utot, elements['ddl'][i],
                                              elements['E'][i], elements['alpha'][i],
                                              elements['dT'][i],
                                              elements['xi'][i], elements['yi'][i],
                                              elements['xj'][i], elements['yj'][i],
                                              elements['yplus'][i], elements['ymoins'][i])))
    tab_sigma_contrainte[i][1] = elements['Sy'][i]/tab_sigma_contrainte[i][0]

print('\nContraintes et F.O.S:\n')
for i in range(nb_element):
    print(f"Élément {i+1}: \u03c3max = {tab_sigma_contrainte[i][0]:.3}{P},\t")
    print(f'FS = %.2f' % tab_sigma_contrainte[i][1])


