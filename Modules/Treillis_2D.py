import numpy as np

from Modules.Fonctions_partagées import (calculer_k_barre2d, assembler_matrice, extraire_matrice, extraire_vecteur,
                                         reconstruire_vecteur, assembler_vecteur)

# Unites : N-mm-MPa
F = 'N'
L = 'mm'
P = 'MPa'


# ---------
# Fonctions
# ---------

# Calcul du vecteur des charges equivalentes dues a la dilatation thermique d'un element Barre2D

def calculer_feq_barre2d(e, a, alpha, dt, xi, yi, xj, yj):
    if dt != 0:
        l_barre = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
        cx = (xj - xi) / l_barre
        cy = (yj - yi) / l_barre
        feq = (e * a * alpha * dt) * np.array([[-cx],
                                               [-cy],
                                               [cx],
                                               [cy]])
    else:
        feq = np.array([[0], [0], [0], [0]])
    return feq


# Calcul de la contrainte dans un element Barre2D
def calculer_contrainte_barre2d(u_tot, ddl, e, alpha, dt, xi, yi, xj, yj):
    ui = u_tot[ddl[0] - 1][0]
    vi = u_tot[ddl[1] - 1][0]
    uj = u_tot[ddl[2] - 1][0]
    vj = u_tot[ddl[3] - 1][0]
    l_ = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / l_
    cy = (yj - yi) / l_
    sig = (e / l_) * np.array([[-cx, -cy, cx, cy]]) @ np.array([[ui], [vi], [uj], [vj]]) - np.array([[e * alpha * dt]])
    return sig[0][0]


# ----------------------------
# Proprietes de chaque element
# ----------------------------


redo = True
while redo is True:
    try:
        nb_noeuds = int(input("Combien de noeuds contient la structure? "))
    except ValueError:
        continue
    noeuds = {'x': [0] * nb_noeuds, 'y': [0] * nb_noeuds, 'ddlx': [0] * nb_noeuds, 'ddly': [0] * nb_noeuds}
    for i in range(nb_noeuds):
        # boucle for pour passer à travers les noeuds, boucle while pour valider les entrées
        while True:
            try:
                noeuds['x'][i] = float(input(f"Position x du noeud {i + 1} en {L}: "))
            except ValueError:
                print('Entrée invalide, corrigez')
                continue
            # commence à 1 quand i est à 0
            noeuds['ddlx'][i] = 2 * i + 1
            try:
                noeuds['y'][i] = float(input(f"Position y du noeud {i + 1} en {L}: "))
            except ValueError or SyntaxError:
                print('Entrée invalide, corrigez')
                continue
            # commence à 2 quand i est à 2
            noeuds['ddly'][i] = 2 * i + 2

            # Si aucun problème d'entrée, sortir de la boucle d'entrée
            break

    print(noeuds)
    redo = bool(input('Appuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

redo = True
while redo is True:
    try:
        nb_elements = int(input("Combien d'éléments contient la structure? "))
    except ValueError or SyntaxError:
        continue
    vide = [0] * int(nb_elements)
    elements = {'ddl': vide.copy(), 'xi': vide.copy(), 'yi': vide.copy(), 'xj': vide.copy(), 'yj': vide.copy(),
                'E': vide.copy(), 'A': vide.copy(), 'alpha': vide.copy(), 'dT': vide.copy(),
                'k': vide.copy(), 'feq': vide.copy()}

    for i in range(nb_elements):
        # soustraction de 1 pour passer du numéro du noeud à son indice dans le tableau
        while True:
            try:
                noeud_i = int(input(f"\n ÉLÉMENT {i + 1}: \nNoeud avant l'élément: ")) - 1
                noeud_j = int(input("Noeud après l'élément: ")) - 1
            except ValueError:
                print("Valeur invalide")
                continue

            elements['ddl'][i] = np.array(
                [noeuds['ddlx'][noeud_i],
                 noeuds['ddly'][noeud_i],
                 noeuds['ddlx'][noeud_j],
                 noeuds['ddly'][noeud_j]])
            elements['xi'][i], elements['yi'][i] = noeuds['x'][noeud_i], noeuds['y'][noeud_i]
            elements['xj'][i], elements['yj'][i] = noeuds['x'][noeud_j], noeuds['y'][noeud_j]

            try:
                print(f"\tPour un ressort, poser un module d'élasticité de 0")
                elements['E'][i] = eval(input(f"Module d'élasticité en {P}: "))
                if elements['E'][i] > 0:
                    elements['A'][i] = eval(input(f"Aire de section en {L}^2: "))
                else:
                    elements['A'][i] = float(input(f"Raideur du ressort en {F}/{L}: "))
                elements['dT'][i] = float(input('Différence de température: '))
                if elements['dT'][i] != 0:
                    elements['alpha'][i] = float(input("Coefficient de dilatation thermique: "))
            except SyntaxError:
                print("Erreur dans les valeurs entrées")
                continue
            except ValueError:
                print("Erreur dans les valeurs entrées")
                continue
            except TypeError:
                print("Erreur dans les valeurs entrées")
                continue
            elements['k'][i] = calculer_k_barre2d(elements['E'][i], elements['A'][i],
                                                  elements['xi'][i], elements['yi'][i],
                                                  elements['xj'][i], elements['yj'][i])
            elements['feq'][i] = calculer_feq_barre2d(elements['E'][i], elements['A'][i],
                                                      elements['alpha'][i], elements['dT'][i],
                                                      elements['xi'][i], elements['yi'][i],
                                                      elements['xj'][i], elements['yj'][i])
            break

    print(elements)
    redo = bool(input('Appuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

Ktot = np.zeros((nb_noeuds * 2, nb_noeuds * 2))
for i in range(nb_elements):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])

Feqtot = np.zeros((nb_noeuds * 2, 1))
for i in range(nb_elements):
    Feqtot = assembler_vecteur(Feqtot, elements['feq'][i], elements['ddl'][i])

# -------------------------
# Conditions aux frontieres
# -------------------------
redo = True
while redo is True:
    nb_Uc = int(input('Combien de déplacements sont connus? '))
    ddlUc = [0] * nb_Uc
    Uc = np.zeros((nb_Uc, 1))
    for i in range(nb_Uc):
        ddlUc[i] = int(input(f'Numéro du ddl connu #{i + 1}: '))
        Uc[i][0] = eval(input(f'Déplacement en {L} du noeud {ddlUc[i]}: '))
    for i in range(nb_Uc):
        print(f'Noeud {ddlUc[i]} : {Uc[i][0]} {L}')
    redo = bool(input('Appuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

redo = True
while redo is True:
    nb_Fc = int(input('Combien de forces sont connues? '))
    ddlFc = [0] * nb_Fc
    Fc = np.zeros((nb_Fc, 1))
    for i in range(nb_Fc):
        ddlFc[i] = int(input(f'Numéro de la force connue #{i + 1}: '))
        Fc[i][0] = eval(input(f'Grandeur en {F} de la force {ddlFc[i]}: '))
    for i in range(nb_Fc):
        print(f'Noeud {ddlFc[i]} : {Fc[i][0]:.2} {F}')
    redo = bool(input('Appuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

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
tab_sigma = [0] * nb_elements
for i in range(nb_elements):
    tab_sigma[i] = calculer_contrainte_barre2d(Utot, elements['ddl'][i], elements['E'][i],
                                               elements['alpha'][i], elements['dT'][i],
                                               elements['xi'][i], elements['yi'][i],
                                               elements['xj'][i], elements['yj'][i])

print(tab_sigma)

tab_force = [0] * nb_elements
for i in range(nb_elements):
    tab_force[i] = elements['A'][i] * tab_sigma[i]

print(tab_force)

for i in range(len(Ui)):
    print(f"Déplacement {ddlFc[i]}:\t{Ui[i]}")
