# Unités : lbf-po
if __name__ == '__main__': pass
import numpy as np

from Modules.Fonctions_partagées import (assembler_matrice, extraire_matrice, extraire_vecteur,
                                         reconstruire_vecteur, noeud_max)


# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Ressort1D
def calculer_k_ressort1d(raideur):
    k = np.array([[raideur, -raideur], [-raideur, raideur]])
    return k


# ----------------------------
# Proprietes de chaque element
# ----------------------------

nb_elements = int(input("Combien d'elements contient la structure? "))
vide = [0]*int(nb_elements)
elements = {'ddl': vide.copy(), 'raideur': vide.copy(), 'k': vide.copy()}

for i in range(nb_elements):
    noeud_1 = input(f"Noeud avant l'élément {i+1}: ")
    noeud_2 = input(f"Noeud après l'élément {i+1}: ")
    elements['ddl'][i] = [int(noeud_1), int(noeud_2)]
    elements['raideur'][i] = int(input(f"Raideur de l'élément {i+1}: "))
    elements['k'][i] = calculer_k_ressort1d(elements['raideur'][i])


# ----------
# Assemblage
# ----------
nb_noeuds = noeud_max(elements['ddl'])

Ktot = np.zeros((nb_noeuds, nb_noeuds))
for i in range(3):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])


# ---------
# Couplages
# ---------

matA = np.array([[0, -2, 0, 1, 0, 0],
                 [0, -3, 0, 0, 0, 1]])

# ------------
# Augmentation
# ------------

Kaug = np.zeros((8, 8))
ddlinit = np.array([1, 2, 3, 4, 5, 6])
ddlaugm = np.array([7, 8])
Kaug = assembler_matrice(Kaug, Ktot, ddlinit, ddlinit)
Kaug = assembler_matrice(Kaug, matA, ddlaugm, ddlinit)
Kaug = assembler_matrice(Kaug, matA.transpose(), ddlinit, ddlaugm)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([1, 3, 5])
Uc = np.array([[0], [0], [0]])
ddlFc = np.array([2, 4, 6, 7, 8])
Fc = np.array([[0], [100], [0], [0], [0]])

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Kaug, ddlFc, ddlFc)
Kcc = extraire_matrice(Kaug, ddlFc, ddlUc)
Kii = extraire_matrice(Kaug, ddlUc, ddlFc)
Kci = extraire_matrice(Kaug, ddlUc, ddlUc)

# --------
# Solution
# --------

Ui = np.linalg.inv(Kic) @ (Fc - Kcc @ Uc)
Fi = Kii @ Ui + Kci @ Uc

# --------------
# Reconstruction
# --------------

# Incluant les multiplicateurs de Lagrange
Utot = reconstruire_vecteur(Uc, ddlUc, Ui, ddlFc)

# Retrait des multiplicateurs de Lagrange
Utot = extraire_vecteur(Utot, ddlinit)

# --------
# Reponses
# --------

all_ressort1 = extraire_vecteur(Utot, [2])
all_ressort2 = extraire_vecteur(Utot, [4])
all_ressort3 = extraire_vecteur(Utot, [6])
print('Allongement du ressort 1 = %.3f po' % all_ressort1[0][0])
print('Allongement du ressort 2 = %.3f po' % all_ressort2[0][0])
print('Allongement du ressort 3 = %.3f po' % all_ressort3[0][0])