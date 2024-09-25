# Unités : lbf-po-psi

import numpy as np

np.set_printoptions(precision=4, linewidth=np.inf)


# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Poutre2D
def calculer_k_poutre2d(E, A, Iz, xi, yi, xj, yj):
    L = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / L
    cy = (yj - yi) / L
    EAL = E * A / L
    EIL = E * Iz / L
    EIL2 = E * Iz / L ** 2
    EIL3 = E * Iz / L ** 3
    klocal = np.array([[EAL, 0, 0, -EAL, 0, 0],
                       [0, 12 * EIL3, 6 * EIL2, 0, -12 * EIL3, 6 * EIL2],
                       [0, 6 * EIL2, 4 * EIL, 0, -6 * EIL2, 2 * EIL],
                       [-EAL, 0, 0, EAL, 0, 0],
                       [0, -12 * EIL3, -6 * EIL2, 0, 12 * EIL3, -6 * EIL2],
                       [0, 6 * EIL2, 2 * EIL, 0, -6 * EIL2, 4 * EIL]])
    R = np.array([[cx, cy, 0, 0, 0, 0],
                  [-cy, cx, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, cx, cy, 0],
                  [0, 0, 0, -cy, cx, 0],
                  [0, 0, 0, 0, 0, 1]])
    kglobal = np.transpose(R) @ klocal @ R
    return kglobal


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


# Reconstruction d'un vecteur a partir de deux vecteurs
def reconstruire_vecteur(Vin1, lig1, Vin2, lig2):
    Vout = np.zeros((len(lig1) + len(lig2), 1))
    for i in range(len(lig1)):
        Vout[lig1[i] - 1][0] = Vin1[i][0]
    for i in range(len(lig2)):
        Vout[lig2[i] - 1][0] = Vin2[i][0]
    return Vout


# Extraction d'un vecteur
def extraire_vecteur(Vin, lig):
    Vout = np.zeros((len(lig), 1))
    for i in range(len(lig)):
        Vout[i][0] = Vin[lig[i] - 1][0]
    return Vout


# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1 = np.array([1, 2, 3, 4, 5, 6])
xi1, yi1 = 0, 0
xj1, yj1 = 12, 0
E1 = 1e7
A1 = 1
Iz1 = 2
k1 = calculer_k_poutre2d(E1, A1, Iz1, xi1, yi1, xj1, yj1)

ddl2 = np.array([4, 5, 6, 7, 8, 9])
xi2, yi2 = 12, 0
xj2, yj2 = 36, 0
E2 = 1e7
A2 = 1
Iz2 = 2
k2 = calculer_k_poutre2d(E2, A2, Iz2, xi2, yi2, xj2, yj2)

ddl3 = np.array([10, 11, 12, 13, 14, 15])
xi3, yi3 = 12, 0
xj3, yj3 = 24, 12
E3 = 1e7
A3 = 0.5
Iz3 = 0.8
k3 = calculer_k_poutre2d(E3, A3, Iz3, xi3, yi3, xj3, yj3)

ddl4 = np.array([13, 14, 15, 16, 17, 18])
xi4, yi4 = 24, 12
xj4, yj4 = 36, 0
E4 = 1e7
A4 = 0.5
Iz4 = 0.8
k4 = calculer_k_poutre2d(E4, A4, Iz4, xi4, yi4, xj4, yj4)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((18, 18))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)
Ktot = assembler_matrice(Ktot, k3, ddl3, ddl3)
Ktot = assembler_matrice(Ktot, k4, ddl4, ddl4)

# ---------
# Couplages
# ---------

matA = np.array([[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0]])

# ------------
# Augmentation
# ------------

Kaug = np.zeros((21, 21))
ddlinit = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18])
ddlaugm = np.array([19, 20, 21])
Kaug = assembler_matrice(Kaug, Ktot, ddlinit, ddlinit)
Kaug = assembler_matrice(Kaug, matA, ddlaugm, ddlinit)
Kaug = assembler_matrice(Kaug, matA.transpose(), ddlinit, ddlaugm)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([1, 2, 3])
Uc = np.array([[0], [0], [0]])
ddlFc = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21])
Fc = np.array([[0], [0], [0], [0], [0], [0], [0], [0], [0], [-500], [0], [0], [0], [0], [0], [0], [0], [0]])

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

# Obtention des forces externes
Ftot = Ktot @ Utot

# --------
# Reponses
# --------

depABC = extraire_vecteur(Utot, [4])
depBDC = extraire_vecteur(Utot, [10])
dep_relatif = abs(depABC[0][0] - depBDC[0][0])
print('Deplacement relatif = %.4f po' % dep_relatif)
