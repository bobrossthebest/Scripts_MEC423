# Unites : N-mm-MPa

import numpy as np

from Modules.Fonctions_partagées import (calculer_k_poutre1d, assembler_matrice, extraire_matrice, extraire_vecteur,
                                         reconstruire_vecteur, assembler_vecteur)


# ---------
# Fonctions
# ---------

# Calcul du vecteur des charges equivalentes dues a une charge repartie constante appliquee sur un element Poutre1D
def calculer_feq_poutre1d(q, L):
    feq = np.array([[q * L / 2],
                    [q * L ** 2 / 12],
                    [q * L / 2],
                    [-q * L ** 2 / 12]])
    return feq


# Calcul par interpolation du deplacement transversal dans un element Poutre1D
def calculer_deplacement_poutre1d(Utot, ddl, L, x):
    vi = Utot[ddl[0] - 1][0]
    ti = Utot[ddl[1] - 1][0]
    vj = Utot[ddl[2] - 1][0]
    tj = Utot[ddl[3] - 1][0]
    dep = (2 * x ** 3 / L ** 3 - 3 * x ** 2 / L ** 2 + 1) * vi + (x ** 3 / L ** 2 - 2 * x ** 2 / L + x) * ti + (
                -2 * x ** 3 / L ** 3 + 3 * x ** 2 / L ** 2) * vj + (x ** 3 / L ** 2 - x ** 2 / L) * tj
    return dep


# Calcul de la contrainte dans un element Poutre1D
def calculer_contrainte_poutre1d(Utot, ddl, E, L, x, y):
    vi = Utot[ddl[0] - 1][0]
    ti = Utot[ddl[1] - 1][0]
    vj = Utot[ddl[2] - 1][0]
    tj = Utot[ddl[3] - 1][0]
    sigma = -y * E * ((12 * x / L ** 3 - 6 / L ** 2) * vi + (6 * x / L ** 2 - 4 / L) * ti + (
                -12 * x / L ** 3 + 6 / L ** 2) * vj + (6 * x / L ** 2 - 2 / L) * tj)
    return sigma


# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1 = np.array([1, 2, 3, 4])
L1 = 400
Iz1 = 384
E1 = 69000
ymax1 = -8
k1 = calculer_k_poutre1d(E1, Iz1, L1)

ddl2 = np.array([3, 4, 5, 6])
L2 = 400
Iz2 = 384
E2 = 69000
q2 = 0.003
ymax2 = -8
k2 = calculer_k_poutre1d(E2, Iz2, L2)
feq2 = calculer_feq_poutre1d(q2, L2)

ddl3 = np.array([5, 6, 7, 8])
L3 = 400
Iz3 = 384
E3 = 69000
ymax3 = -8
k3 = calculer_k_poutre1d(E3, Iz3, L3)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((8, 8))
Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)
Ktot = assembler_matrice(Ktot, k3, ddl3, ddl3)

Feqtot = np.zeros((8, 1))
Feqtot = assembler_vecteur(Feqtot, feq2, ddl2)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([1, 7, 8])
Uc = np.array([[0], [0], [0]])
ddlFc = np.array([2, 3, 4, 5, 6])
Fc = np.array([[0], [0], [0], [0], [0]])

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

xB = 200
V_B = calculer_deplacement_poutre1d(Utot, ddl1, L1, xB)
print('V_B = %.3f mm' % V_B)

R_A = extraire_vecteur(Ftot, [1])
print('R_A = %.3f N' % R_A[0][0])

xE = 400
sigmaMax_E = calculer_contrainte_poutre1d(Utot, ddl3, E3, L3, xE, ymax3)
print('SigmaMax_E = %.2f MPa' % sigmaMax_E)