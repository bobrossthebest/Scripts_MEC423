# Unités : N-mm-MPa

import numpy as np

from Modules.Fonctions_partagées import (calculer_k_poutre1d, assembler_matrice, extraire_matrice)


# ---------
# Fonctions
# ---------


# Calcul de la matrice de rigidite non-lineaire d'un element Poutre1D
def calculer_knl_poutre1d(L):
    knl = np.array([[1.2 / L, 0.1, -1.2 / L, 0.1],
                    [0.1, 2 * L / 15, -0.1, -L / 30],
                    [-1.2 / L, -0.1, 1.2 / L, -0.1],
                    [0.1, -L / 30, -0.1, 2 * L / 15]])
    return knl


# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1 = np.array([1, 2, 3, 4])
L1 = 250
Iz1 = 12 * 12 ** 3 / 12
E1 = 2e5
k1 = calculer_k_poutre1d(E1, Iz1, L1)
knl1 = calculer_knl_poutre1d(L1)

ddl2 = np.array([3, 4, 5, 6])
L2 = 200
Iz2 = np.pi * 20 ** 4 / 64
E2 = 2e5
k2 = calculer_k_poutre1d(E2, Iz2, L2)
knl2 = calculer_knl_poutre1d(L2)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((6, 6))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)

Knltot = np.zeros((6, 6))
Knltot = assembler_matrice(Knltot, knl1, ddl1, ddl1)
Knltot = assembler_matrice(Knltot, knl2, ddl2, ddl2)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([5, 6])
Uc = np.array([[0], [0]])
ddlFc = np.array([1, 2, 3, 4])
Fc = np.array([[0], [0], [0], [0]])

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Ktot, ddlFc, ddlFc)
Knlic = extraire_matrice(Knltot, ddlFc, ddlFc)

# --------------------------------------
# Solution (valeurs et vecteurs propres)
# --------------------------------------

valprop, vecprop = np.linalg.eig(-np.linalg.inv(Knlic) @ Kic)

# --------
# Reponses
# --------

# Identification de la charge critique
col = -1
Pcr = -np.inf
for i in range(len(Fc)):
    if valprop[i] < 0 and valprop[i] > Pcr:
        col = i
        Pcr = valprop[i]
print('Charge critique = %.0f N (colonne %d)' % (Pcr, col + 1))

# Premier mode de flambage
Ui = vecprop[:, col]
print('Valeur des DDL libres du 1er mode de flambage :')
print(Ui)
