# ---------
# imports externes
# ---------

import numpy as np


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


# Assemblage d'un vecteur
def assembler_vecteur(Ftot, f, lig):
    for i in range(len(lig)):
        Ftot[lig[i] - 1][0] += f[i][0]
    return Ftot


# Calcul de la matrice de rigidite d'un element Poutre1D
def calculer_k_poutre1d(E, I, L):
    k = (E * I / L ** 3) * np.array([[12, 6 * L, -12, 6 * L],
                                     [6 * L, 4 * L ** 2, -6 * L, 2 * L ** 2],
                                     [-12, -6 * L, 12, -6 * L],
                                     [6 * L, 2 * L ** 2, -6 * L, 4 * L ** 2]])
    return k


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
