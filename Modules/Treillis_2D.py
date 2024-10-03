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
nb_noeuds = int(input("Combien de noeuds contient la structure? "))
noeuds = {'x': [0]*nb_noeuds, 'y': [0]*nb_noeuds, 'ddlx': [0]*nb_noeuds, 'ddly': [0]*nb_noeuds}
for i in range(nb_noeuds):
    # Le premier noeud, d'indice 0, est le noeud "1" à l'affichage
    noeuds['x'][i] = float(input(f"Position x du noeud {i+1} en {F}: "))
    # commence à 1 quand i est à 0
    noeuds['ddlx'][i] = 2*i+1
    noeuds['y'][i] = float(input(f"Position y du noeud {i+1} en {F}: "))
    # commence à 2 quand i est à 2
    noeuds['ddly'][i] = 2*i+2

nb_elements = int(input("\nCombien d'éléments contient la structure? "))
vide = [0]*int(nb_elements)
elements = {'ddl': vide.copy(), 'xi': vide.copy(), 'yi': vide.copy(), 'xj': vide.copy(), 'yj': vide.copy(),
            'E': vide.copy(), 'A': vide.copy(), 'alpha': vide.copy(), 'dT': vide.copy(),
            'k': vide.copy(), 'feq': vide.copy()}

for i in range(nb_elements):
    # soustraction de 1 pour passer du numéro du noeud à son indice dans le tableau
    noeud_i = int(input(f"\n ÉLÉMENT {i+1}: \nNoeud avant l'élément: ")) - 1
    noeud_j = int(input("Noeud après l'élément: ")) - 1
    elements['ddl'][i] = np.array(
        [noeuds['ddlx'][noeud_i],
         noeuds['ddly'][noeud_i],
         noeuds['ddlx'][noeud_j],
         noeuds['ddlx'][noeud_j]])
    elements['xi'][i], elements['yi'][i] = noeuds['x'][noeud_i], noeuds['y'][noeud_i]
    elements['xj'][i], elements['yj'][i] = noeuds['x'][noeud_j], noeuds['y'][noeud_j]

    print(f"\tPour un ressort, poser un module d'élasticité de 0 et la raideur au-lieu de l'aire")
    elements['E'][i] = float(input(f"Module d'élasticité en {P}: "))
    if elements['E'][i] > 0:
        elements['A'][i] = float(input(f"Aire de section en {L}^2: "))
    elements['dT'][i] = float(input('Différence de température: '))
    if elements['dT'][i] != 0:
        elements['alpha'][i] = float(input("Coefficient de dilatation thermique: "))
    elements['k'][i] = calculer_k_barre2d(elements['E'][i], elements['A'][i],
                                          elements['xi'][i], elements['yi'][i],
                                          elements['xj'][i], elements['yj'][i])
    elements['feq'][i] = calculer_feq_barre2d(elements['E'][i], elements['A'][i],
                                              elements['alpha'][i], elements['dT'][i],
                                              elements['xi'][i], elements['yi'][i],
                                              elements['xj'][i], elements['yj'][i])


Ktot = np.zeros((nb_noeuds * 2, nb_noeuds * 2))
for i in range(nb_elements):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])

Feqtot = np.zeros((nb_noeuds * 2, 1))
for i in range(nb_elements):
    Feqtot = assembler_vecteur(Feqtot, elements['feq'][i], elements['ddl'][i])

# -------------------------
# Conditions aux frontieres
# -------------------------
nb_Uc = int(input('Combien de déplacements sont connus?'))
ddlUc = [0]*nb_Uc
Uc = np.zeros((nb_Uc, 1))
for i in range(nb_Uc):
    ddlUc[i] = int(input(f'Numéro du ddl connu #{i+1}: '))
    Uc[i][0] = eval(input(f'Déplacement en {L} du noeud {ddlUc[i]}: '))
nb_Fc = int(input('Combien de forces sont connues?'))
ddlFc = [0]*nb_Fc
Fc = np.zeros((nb_Fc, 1))
for i in range(nb_Fc):
    ddlFc[i] = int(input(f'Numéro de la force connue #{i+1}: '))
    Fc[i][0] = eval(input(f'Grandeur en {F} de la force {ddlFc[i]}: '))

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

VF = extraire_vecteur(Utot, [8])
print('UyF = %.2f mm' % VF[0][0])

RB = extraire_vecteur(Ftot, [1, 2])
print('RxB = %.1f N, RyB = %.1f N' % (RB[0][0], RB[1][0]))

tab_sigma = [0] * nb_elements
for i in range(nb_elements):
    tab_sigma[i] = calculer_contrainte_barre2d(Utot, elements['ddl'][i], elements['E'][i],
                                               elements['alpha'][i], elements['dT'][i],
                                               elements['xi'][i], elements['yi'][i],
                                               elements['xj'][i], elements['yj'][i])

tab_force = [0]*nb_elements
for i in range(nb_elements):
    tab_force[i] = elements['A'][i]*tab_sigma[i]

print(f'force_CF = {tab_force[2]:.1f} {F}')
