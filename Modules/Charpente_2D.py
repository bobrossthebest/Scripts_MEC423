# Unites : N-mm-MPa

import numpy as np
from Modules.Fonctions_partagées import (assembler_matrice, extraire_matrice, extraire_vecteur, reconstruire_vecteur,
                                         calculer_k_poutre2d, calcul_Iz)

# Unites :
F = input("\nQuelle est l'unité de mesure de force?\t\t")
L = input("Quelle est l'unité de mesure de longueur?\t")
P = input("Quelle est l'unité de mesure de contrainte?\t")

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
                'xi': vide.copy(),
                'yi': vide.copy(),
                'xj': vide.copy(),
                'yj': vide.copy(),
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


            try:
                elements['Sy'][i] = eval(input(f"Limite d'ecoulement en {P}:\t"))
                elements['Iz'][i] = calcul_Iz(L)
                elements['yplus'][i] = eval(input(f"Y plus de la section en {L}:\t"))
                elements['ymoins'][i] = eval(input(f"Y moins de la section en {L}:\t"))
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

ddl1 = np.array([1, 2, 3, 4, 5, 6])
xi1, yi1 = 0, 0
xj1, yj1 = 4000, 3000
A1 = 24 * 24
Iz1 = 24 * 24 ** 3 / 12
E1 = 2e5
Sy1 = 250
yplus1 = 12
ymoins1 = -12
qx1 = 0.15
qy1 = -0.2
alpha1 = 12e-6
dT1 = 0
k1 = calculer_k_poutre2d(E1, A1, Iz1, xi1, yi1, xj1, yj1)
feq1 = calculer_feq_poutre2d(E1, A1, alpha1, dT1, qx1, qy1, xi1, yi1, xj1, yj1)

ddl2 = np.array([4, 5, 6, 7, 8, 9])
xi2, yi2 = 4000, 3000
xj2, yj2 = 4000, 0
A2 = 24 * 24
Iz2 = 24 * 24 ** 3 / 12
E2 = 2e5
Sy2 = 250
yplus2 = 12
ymoins2 = -12
qx2 = 0
qy2 = 0
alpha2 = 12e-6
dT2 = 40
k2 = calculer_k_poutre2d(E2, A2, Iz2, xi2, yi2, xj2, yj2)
feq2 = calculer_feq_poutre2d(E2, A2, alpha2, dT2, qx2, qy2, xi2, yi2, xj2, yj2)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((9, 9))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)

Feqtot = np.zeros((9, 1))
Feqtot = assembler_vecteur(Feqtot, feq1, ddl1)
Feqtot = assembler_vecteur(Feqtot, feq2, ddl2)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([1, 2, 3, 7, 8])
Uc = np.array([[0], [0], [0], [0], [0]])
ddlFc = np.array([4, 5, 6, 9])
Fc = np.array([[0], [0], [0], [0]])

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Ktot, ddlFc, ddlFc)
Kcc = extraire_matrice(Ktot, ddlFc, ddlUc)
Kii = extraire_matrice(Ktot, ddlUc, ddlFc)
Kci = extraire_matrice(Ktot, ddlUc, ddlUc)
Feqi = extraire_vecteur(Feqtot, ddlUc)
Feqc = extraire_vecteur(Feqtot, ddlFc)

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

UxB = extraire_vecteur(Utot, [4])
print('UxB = %.3f mm' % UxB[0][0])

sig1 = calculer_contraintes_poutre2d(Utot, ddl1, E1, alpha1, dT1, xi1, yi1, xj1, yj1, yplus1, ymoins1)
sigMax1 = max(abs(sig1))
FS1 = Sy1 / sigMax1
print('FS1 = %.2f' % FS1)
