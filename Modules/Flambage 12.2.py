# Unités : N-mm-MPa

import numpy as np


# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Poutre1D
def calculer_k_poutre1d(E, I, L):
    k = (E * I / L ** 3) * np.array([[12, 6 * L, -12, 6 * L],
                                     [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
                                     [-12, -6 * L, 12, -6 * L],
                                     [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])
    return k


# Calcul de la matrice de rigidite non-lineaire d'un element Poutre1D
def calculer_knl_poutre1d(L):
    knl = np.array([[1.2 / L, 0.1, -1.2 / L, 0.1],
                    [0.1, 2 * L / 15, -0.1, -L / 30],
                    [-1.2 / L, -0.1, 1.2 / L, -0.1],
                    [0.1, -L / 30, -0.1, 2 * L / 15]])
    return knl


# Assemblage d'une matrice
def assembler_matrice(Ktot, k, lig, col):
    for i in range(len(lig)):
        for j in range(len(col)):
            Ktot[lig[i] - 1][col[j] - 1] += k[i][j]
    return Ktot


# Extraction d'une sous-matrice
def extraire_matrice(Kin, lig, col):
    Kout = np.zeros((len(lig), len(col)))
    for i in range(len(lig)):
        for j in range(len(col)):
            Kout[i][j] = Kin[lig[i] - 1][col[j] - 1]
    return Kout


# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1 = np.array([1, 2, 3, 4])
L1 = 300
Iz1 = 15**4/12
E1 = 69e3
k1 = calculer_k_poutre1d(E1, Iz1, L1)
knl1 = calculer_knl_poutre1d(L1)

ddl2 = np.array([3, 4, 5, 6])
L2 = 300
Iz2 = Iz1
E2 = E1
k2 = calculer_k_poutre1d(E2, Iz2, L2)
knl2 = calculer_knl_poutre1d(L2)

ddl3 = np.array([5, 6, 7, 8])
L3 = 400
Iz3 = 20**4/12
E3 = 210e3
k3 = calculer_k_poutre1d(E3, Iz3, L3)
knl3 = calculer_knl_poutre1d(L3)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((8, 8))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)
Ktot = assembler_matrice(Ktot, k3, ddl3, ddl3)

Knltot = np.zeros((8, 8))
Knltot = assembler_matrice(Knltot, knl1, ddl1, ddl1)
Knltot = assembler_matrice(Knltot, knl2, ddl2, ddl2)
Knltot = assembler_matrice(Knltot, knl3, ddl3, ddl3)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([3, 7])
Uc = np.array([[0], [0]])
ddlFc = np.array([1, 2, 4, 5, 6, 8])
Fc = np.array([[0], [0], [0], [0], [0], [0]])

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
print('\nCharge critique = %.0f N (colonne %d)\n' % (Pcr, col + 1))

# Premier mode de flambage
Ui = vecprop[:, col]
print('Valeur des DDL libres du 1er mode de flambage :')
print(Ui)

P = -1200

print(f"\nFacteur de sécurité p/rapport à P : {Pcr/P:.3}")
