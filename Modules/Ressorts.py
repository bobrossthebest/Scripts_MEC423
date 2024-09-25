# Unités : lbf-po

import numpy as np

from Modules.Fonctions_partagées import (assembler_matrice, extraire_matrice, extraire_vecteur,
                                         reconstruire_vecteur)


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

ddl1 = np.array([1, 2])
raid1 = 200
k1 = calculer_k_ressort1d(raid1)

ddl2 = np.array([3, 4])
raid2 = 120
k2 = calculer_k_ressort1d(raid2)

ddl3 = np.array([5, 6])
raid3 = 80
k3 = calculer_k_ressort1d(raid3)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((6, 6))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)
Ktot = assembler_matrice(Ktot, k3, ddl3, ddl3)

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
