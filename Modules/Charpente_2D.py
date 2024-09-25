# Unites : N-mm-MPa

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


# Calcul du vecteur des charges equivalentes dues a une charge repartie constante appliquee sur un element Poutre2D
def calculer_feq_poutre2d(E, A, alpha, dT, qx, qy, xi, yi, xj, yj):
    L = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / L
    cy = (yj - yi) / L
    feqT = A * E * alpha * dT * np.array([[-cx],
                                          [-cy],
                                          [0],
                                          [cx],
                                          [cy],
                                          [0]])
    feqQ = np.array([[qx * L / 2],
                     [qy * L / 2],
                     [(qy * cx - qx * cy) * L ** 2 / 12],
                     [qx * L / 2],
                     [qy * L / 2],
                     [(qx * cy - qy * cx) * L ** 2 / 12]])
    feq = feqT + feqQ
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


# Calcul des contraintes aux extremites d'un element Poutre2D
def calculer_contraintes_poutre2d(Utot, ddl, E, alpha, dT, xi, yi, xj, yj, yplus, ymoins):
    L = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / L
    cy = (yj - yi) / L
    R = np.array([[cx, cy, 0, 0, 0, 0],
                  [-cy, cx, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, cx, cy, 0],
                  [0, 0, 0, -cy, cx, 0],
                  [0, 0, 0, 0, 0, 1]])
    dispGlobal = np.array([[Utot[ddl[0] - 1][0]],
                           [Utot[ddl[1] - 1][0]],
                           [Utot[ddl[2] - 1][0]],
                           [Utot[ddl[3] - 1][0]],
                           [Utot[ddl[4] - 1][0]],
                           [Utot[ddl[5] - 1][0]]])
    dispLocal = R @ dispGlobal
    ui = dispLocal[0][0]
    vi = dispLocal[1][0]
    ti = dispLocal[2][0]
    uj = dispLocal[3][0]
    vj = dispLocal[4][0]
    tj = dispLocal[5][0]
    sigTrac = E * (uj - ui) / L - E * alpha * dT
    sigFlexiPlus = -yplus * E * (6 * (vj - vi) / L ** 2 - 2 * (2 * ti + tj) / L)
    sigFlexiMoins = -ymoins * E * (6 * (vj - vi) / L ** 2 - 2 * (2 * ti + tj) / L)
    sigFlexjPlus = -yplus * E * (6 * (vi - vj) / L ** 2 + 2 * (ti + 2 * tj) / L)
    sigFlexjMoins = -ymoins * E * (6 * (vi - vj) / L ** 2 + 2 * (ti + 2 * tj) / L)
    sigCombiPlus = sigFlexiPlus + sigTrac
    sigCombiMoins = sigFlexiMoins + sigTrac
    sigCombjPlus = sigFlexjPlus + sigTrac
    sigCombjMoins = sigFlexjMoins + sigTrac
    sig = np.array(
        [sigTrac, sigFlexiPlus, sigFlexiMoins, sigFlexjPlus, sigFlexjMoins, sigCombiPlus, sigCombiMoins, sigCombjPlus,
         sigCombjMoins])
    return sig


# ----------------------------
# Proprietes de chaque element
# ----------------------------

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
