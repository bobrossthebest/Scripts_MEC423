# Unites : N-mm-MPa

import numpy as np

np.set_printoptions(precision=4, linewidth=np.inf)


# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Barre2D


def calculer_k_barre2d(E, A, xi, yi, xj, yj):
    L = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / L
    cy = (yj - yi) / L
    k = (E * A / L) * np.array([[cx ** 2, cx * cy, -cx ** 2, -cx * cy],
                                [cx * cy, cy ** 2, -cx * cy, -cy ** 2],
                                [-cx ** 2, -cx * cy, cx ** 2, cx * cy],
                                [-cx * cy, -cy ** 2, cx * cy, cy ** 2]])
    return k


# Calcul du vecteur des charges equivalentes dues a la dilatation thermique d'un element Barre2D
def calculer_feq_barre2d(E, A, alpha, dT, xi, yi, xj, yj):
    L = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / L
    cy = (yj - yi) / L
    feq = (E * A * alpha * dT) * np.array([[-cx],
                                           [-cy],
                                           [cx],
                                           [cy]])
    return feq


# Assemblage d'une matrice
def assembler_matrice(Ktot, k, lig, col):
    for i in range(len(lig)):
        for j in range(len(col)):
            Ktot[lig[i] - 1][col[j] - 1] += k[i][j]
    return Ktot


# Assemblage d'un vecteur
def assembler_vecteur(Ftot, f, lig):
    for i in range(len(lig)):
        Ftot[lig[i] - 1][0] += f[i][0]
    return Ftot


# Extraction d'une sous-matrice
def extraire_matrice(Kin, lig, col):
    Kout = np.zeros((len(lig), len(col)))
    for i in range(len(lig)):
        for j in range(len(col)):
            Kout[i][j] = Kin[lig[i] - 1][col[j] - 1]
    return Kout


# Extraction d'un vecteur
def extraire_vecteur(Vin, lig):
    Vout = np.zeros((len(lig), 1))
    for i in range(len(lig)):
        Vout[i][0] = Vin[lig[i] - 1][0]
    return Vout


# Reconstruction d'un vecteur a partir de deux vecteurs
def reconstruire_vecteur(Vin1, lig1, Vin2, lig2):
    Vout = np.zeros((len(lig1) + len(lig2), 1))
    for i in range(len(lig1)):
        Vout[lig1[i] - 1][0] = Vin1[i][0]
    for i in range(len(lig2)):
        Vout[lig2[i] - 1][0] = Vin2[i][0]
    return Vout


# Calcul de la contrainte dans un element Barre2D
def calculer_contrainte_barre2d(Utot, ddl, E, alpha, dT, xi, yi, xj, yj):
    Ui = Utot[ddl[0] - 1][0]
    Vi = Utot[ddl[1] - 1][0]
    Uj = Utot[ddl[2] - 1][0]
    Vj = Utot[ddl[3] - 1][0]
    L = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / L
    cy = (yj - yi) / L
    sigma = (E / L) * np.array([[-cx, -cy, cx, cy]]) @ np.array([[Ui], [Vi], [Uj], [Vj]]) - np.array([[E * alpha * dT]])
    return sigma[0][0]


# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1 = np.array([1, 2, 3, 4])
xi1, yi1 = 0, 0
xj1, yj1 = 800, 600
E1 = 2e5
A1 = 24
alpha1 = 12e-6
dT1 = 100
Sy1 = 310
k1 = calculer_k_barre2d(E1, A1, xi1, yi1, xj1, yj1)
feq1 = calculer_feq_barre2d(E1, A1, alpha1, dT1, xi1, yi1, xj1, yj1)

ddl2 = np.array([3, 4, 5, 6])
xi2, yi2 = 800, 600
xj2, yj2 = 800, 0
E2 = 2e5
A2 = 24
alpha2 = 12e-6
dT2 = 0
Sy2 = 310
k2 = calculer_k_barre2d(E2, A2, xi2, yi2, xj2, yj2)
feq2 = calculer_feq_barre2d(E2, A2, alpha2, dT2, xi2, yi2, xj2, yj2)

ddl3 = np.array([3, 4, 7, 8])
xi3, yi3 = 800, 600
xj3, yj3 = 1600, 0
E3 = 2e5
A3 = 24
alpha3 = 12e-6
dT3 = 0
Sy3 = 310
k3 = calculer_k_barre2d(E3, A3, xi3, yi3, xj3, yj3)
feq3 = calculer_feq_barre2d(E3, A3, alpha3, dT3, xi3, yi3, xj3, yj3)

ddl4 = np.array([1, 2, 5, 6])
xi4, yi4 = 0, 0
xj4, yj4 = 800, 0
E4 = 2e5
A4 = 24
alpha4 = 12e-6
dT4 = 0
Sy4 = 310
k4 = calculer_k_barre2d(E4, A4, xi4, yi4, xj4, yj4)
feq4 = calculer_feq_barre2d(E4, A4, alpha4, dT4, xi4, yi4, xj4, yj4)

ddl5 = np.array([5, 6, 7, 8])
xi5, yi5 = 800, 0
xj5, yj5 = 1600, 0
E5 = 2e5
A5 = 24
alpha5 = 12e-6
dT5 = 0
Sy5 = 310
k5 = calculer_k_barre2d(E5, A5, xi5, yi5, xj5, yj5)
feq5 = calculer_feq_barre2d(E5, A5, alpha5, dT5, xi5, yi5, xj5, yj5)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((8, 8))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)
Ktot = assembler_matrice(Ktot, k3, ddl3, ddl3)
Ktot = assembler_matrice(Ktot, k4, ddl4, ddl4)
Ktot = assembler_matrice(Ktot, k5, ddl5, ddl5)

Feqtot = np.zeros((8, 1))
Feqtot = assembler_vecteur(Feqtot, feq1, ddl1)
Feqtot = assembler_vecteur(Feqtot, feq2, ddl2)
Feqtot = assembler_vecteur(Feqtot, feq3, ddl3)
Feqtot = assembler_vecteur(Feqtot, feq4, ddl4)
Feqtot = assembler_vecteur(Feqtot, feq5, ddl5)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([1, 2, 6])
Uc = np.array([[0], [0], [0]])
ddlFc = np.array([3, 4, 5, 7, 8])
Fc = np.array([[0], [0], [0], [-1200 * np.sin(40 * np.pi / 180)], [-1200 * np.cos(40 * np.pi / 180)]])

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

sigma3 = calculer_contrainte_barre2d(Utot, ddl3, E3, alpha3, dT3, xi3, yi3, xj3, yj3)
force3 = sigma3 * A3
print('force_CF = %.1f N' % force3)

sigma1 = calculer_contrainte_barre2d(Utot, ddl1, E1, alpha1, dT1, xi1, yi1, xj1, yj1)
fs1 = Sy1 / sigma1
print('FS_BC = %.2f' % fs1)
