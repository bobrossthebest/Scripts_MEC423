# Unites : lb-po-psi

# Pour des situations où les déplacements sont déjà donnés

import numpy as np

# ---------
# Fonctions
# ---------

# Calcul de la matrice ksi en etat plan de contrainte
def calculer_ksi_EPC(E, nu):
    ksi = E / (1 - nu ** 2) * np.array([[1, nu, 0],
                                        [nu, 1, 0],
                                        [0, 0, (1 - nu) / 2]])
    return ksi


# Calcul de l'aire d'un element T3
def calculer_aire_T3(xi, yi, xj, yj, xk, yk):
    temp = np.array([[1, xi, yi],
                     [1, xj, yj],
                     [1, xk, yk]])
    aire = 0.5 * np.linalg.det(temp)
    return aire


# Calcul de la matrice de B d'un element T3
def calculer_B_T3(xi, yi, xj, yj, xk, yk):
    aire = calculer_aire_T3(xi, yi, xj, yj, xk, yk)
    B = 1 / (2 * aire) * np.array([[yj - yk, 0, yk - yi, 0, yi - yj, 0],
                                   [0, xk - xj, 0, xi - xk, 0, xj - xi],
                                   [xk - xj, yj - yk, xi - xk, yk - yi, xj - xi, yi - yj]])
    return B


# Calcul de la matrice de rigidite d'un element T3
def calculer_k_T3(E, nu, t, xi, yi, xj, yj, xk, yk):
    ksi = calculer_ksi_EPC(E, nu)
    A = calculer_aire_T3(xi, yi, xj, yj, xk, yk)
    B = calculer_B_T3(xi, yi, xj, yj, xk, yk)
    k = t * A * (np.transpose(B) @ ksi @ B)
    return k


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


# Evaluation de la fonction d'interpolation Ni de l'element T3 a un point donne
def evaluer_Ni_T3(xi, yi, xj, yj, xk, yk, x, y):
    aire = calculer_aire_T3(xi, yi, xj, yj, xk, yk)
    Ni = 1 / (2 * aire) * ((xj * yk - xk * yj) + (yj - yk) * x + (xk - xj) * y)
    return Ni


# Evaluation de la fonction d'interpolation Nj de l'element T3 a un point donne
def evaluer_Nj_T3(xi, yi, xj, yj, xk, yk, x, y):
    aire = calculer_aire_T3(xi, yi, xj, yj, xk, yk)
    Nj = 1 / (2 * aire) * ((xk * yi - xi * yk) + (yk - yi) * x + (xi - xk) * y)
    return Nj


# Evaluation de la fonction d'interpolation Nk de l'element T3 a un point donne
def evaluer_Nk_T3(xi, yi, xj, yj, xk, yk, x, y):
    aire = calculer_aire_T3(xi, yi, xj, yj, xk, yk)
    Nk = 1 / (2 * aire) * ((xi * yj - xj * yi) + (yi - yj) * x + (xj - xi) * y)
    return Nk


# Extraction d'un vecteur
def extraire_vecteur(Vin, lig):
    Vout = np.zeros((len(lig), 1))
    for i in range(len(lig)):
        Vout[i][0] = Vin[lig[i] - 1][0]
    return Vout


# Calcul des contraintes dans un element T3
def calculer_contraintes_T3(Utot, ddl, E, nu, xi, yi, xj, yj, xk, yk):
    ksi = calculer_ksi_EPC(E, nu)
    B = calculer_B_T3(xi, yi, xj, yj, xk, yk)
    U = extraire_vecteur(Utot, ddl)
    sig = ksi @ B @ U
    return sig


# Calcul de la contrainte equivalente de Mises en etat plan de contrainte
def calculer_mises_EPC(sig):
    sx = sig[0][0]
    sy = sig[1][0]
    txy = sig[2][0]
    mises = (sx ** 2 + sy ** 2 + 3 * txy ** 2 - sx * sy) ** 0.5
    return mises


# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1 = np.array([1, 2, 3, 4, 5, 6])
xi1, yi1 = 0, 0.787
xj1, yj1 = 0.301, 0.727
xk1, yk1 = 0, 1.083
t1 = 0.197
E1 = 3.9 * 10**6
nu1 = 0.33
k1 = calculer_k_T3(E1, nu1, t1, xi1, yi1, xj1, yj1, xk1, yk1)

# ddl2 = np.array([3, 4, 5, 6, 7, 8])
# xi2, yi2 = 0, 10
# xj2, yj2 = 20, 0
# xk2, yk2 = 20, 10
# t2 = 1.2
# E2 = 2e5
# nu2 = 0.3
# k2 = calculer_k_T3(E2, nu2, t2, xi2, yi2, xj2, yj2, xk2, yk2)

# ----------
# Assemblage
# ----------

# Ktot = np.zeros((8, 8))
# Ktot = assembler_matrice(Ktot, k1, ddl1, ddl1)
# Ktot = assembler_matrice(Ktot, k2, ddl2, ddl2)

# -------------------------
# Conditions aux frontieres
# -------------------------

# ddlUc = np.array([1, 2, 3, 4])
# Uc = np.array([[0], [0], [0], [0]])
# ddlFc = np.array([5, 6, 7, 8])
# Fc = np.array([[0], [0], [0], [-3000]])

# ---------------
# Partitionnement
# ---------------
#
# Kic = extraire_matrice(Ktot, ddlFc, ddlFc)
# Kcc = extraire_matrice(Ktot, ddlFc, ddlUc)
# Kii = extraire_matrice(Ktot, ddlUc, ddlFc)
# Kci = extraire_matrice(Ktot, ddlUc, ddlUc)

# --------
# Solution
# --------

# Ui = np.linalg.inv(Kic) @ (Fc - Kcc @ Uc)
# Fi = Kii @ Ui + Kci @ Uc

# --------------
# Reconstruction
# --------------
# Pour l'exo 6.1
Utot = [[0], [-0.00251], [0.00198],[-0.00228], [0], [-0.00318]]
# Utot = reconstruire_vecteur(Uc, ddlUc, Ui, ddlFc)
# Ftot = reconstruire_vecteur(Fc, ddlFc, Fi, ddlUc)

# --------
# Reponses
# --------

# xP = 17
# yP = 8
# Ni2P = evaluer_Ni_T3(xi2, yi2, xj2, yj2, xk2, yk2, xP, yP)
# Nj2P = evaluer_Nj_T3(xi2, yi2, xj2, yj2, xk2, yk2, xP, yP)
# Nk2P = evaluer_Nk_T3(xi2, yi2, xj2, yj2, xk2, yk2, xP, yP)
# vi2 = extraire_vecteur(Utot, [4])
# vj2 = extraire_vecteur(Utot, [6])
# vk2 = extraire_vecteur(Utot, [8])
# vP = vi2[0][0] * Ni2P + vj2[0][0] * Nj2P + vk2[0][0] * Nk2P
# print('vP = %.4f mm' % vP)

sig1 = calculer_contraintes_T3(Utot, ddl1, E1, nu1, xi1, yi1, xj1, yj1, xk1, yk1)
mises1 = calculer_mises_EPC(sig1)
print('Mises_1 = %.1f MPa' % mises1)
print('sig_x = %.1f MPa' % sig1[0][0])

# sig2 = calculer_contraintes_T3(Utot, ddl2, E2, nu2, xi2, yi2, xj2, yj2, xk2, yk2)
# mises2 = calculer_mises_EPC(sig2)
# print('Mises_2 = %.1f MPa' % mises2)
