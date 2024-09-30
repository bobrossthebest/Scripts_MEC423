# ---------
# imports externes
# ---------

import numpy as np


# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Barre2D
def calculer_k_barre2d(e, a, xi, yi, xj, yj):
    l_barre = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / l_barre
    cy = (yj - yi) / l_barre
    tab_cos =  np.array([[cx ** 2, cx * cy, -cx ** 2, -cx * cy],
                                      [cx * cy, cy ** 2, -cx * cy, -cy ** 2],
                                      [-cx ** 2, -cx * cy, cx ** 2, cx * cy],
                                      [-cx * cy, -cy ** 2, cx * cy, cy ** 2]])
    if e > 0:
        k = (e * a / l_barre) * tab_cos
    else:
        k = a * tab_cos
    return k


# Assemblage d'une matrice
def assembler_matrice(k_tot, k, lig, col):
    for i in range(len(lig)):
        for j in range(len(col)):
            k_tot[lig[i] - 1][col[j] - 1] += k[i][j]
    return k_tot


# Extraction d'une sous-matrice
def extraire_matrice(k_in, lig, col):
    k_out = np.zeros((len(lig), len(col)))
    for i in range(len(lig)):
        for j in range(len(col)):
            k_out[i][j] = k_in[lig[i] - 1][col[j] - 1]
    return k_out


# Reconstruction d'un vecteur a partir de deux vecteurs
def reconstruire_vecteur(v_in1, lig1, v_in2, lig2):
    v_out = np.zeros((len(lig1) + len(lig2), 1))
    for i in range(len(lig1)):
        v_out[lig1[i] - 1][0] = v_in1[i][0]
    for i in range(len(lig2)):
        v_out[lig2[i] - 1][0] = v_in2[i][0]
    return v_out


# Extraction d'un vecteur
def extraire_vecteur(v_in, lig):
    v_out = np.zeros((len(lig), 1))
    for i in range(len(lig)):
        v_out[i][0] = v_in[lig[i] - 1][0]
    return v_out


# Assemblage d'un vecteur
def assembler_vecteur(f_tot, f, lig):
    for i in range(len(lig)):
        f_tot[lig[i] - 1][0] += f[i][0]
    return f_tot


# Calcul de la matrice de rigidite d'un element Poutre1D
def calculer_k_poutre1d(e, i, l_poutre):
    k = (e * i / l_poutre ** 3) * np.array([[12, 6 * l_poutre, -12, 6 * l_poutre],
                                            [6 * l_poutre, 4 * l_poutre ** 2, -6 * l_poutre, 2 * l_poutre ** 2],
                                            [-12, -6 * l_poutre, 12, -6 * l_poutre],
                                            [6 * l_poutre, 2 * l_poutre ** 2, -6 * l_poutre, 4 * l_poutre ** 2]])
    return k


# Calcul de la matrice de rigidite d'un element Poutre2D
def calculer_k_poutre2d(e, a, iz, xi, yi, xj, yj):
    l_poutre = ((xj - xi) ** 2 + (yj - yi) ** 2) ** 0.5
    cx = (xj - xi) / l_poutre
    cy = (yj - yi) / l_poutre
    eal = e * a / l_poutre
    eil = e * iz / l_poutre
    eil2 = e * iz / l_poutre ** 2
    eil3 = e * iz / l_poutre ** 3
    klocal = np.array([[eal, 0, 0, -eal, 0, 0],
                       [0, 12 * eil3, 6 * eil2, 0, -12 * eil3, 6 * eil2],
                       [0, 6 * eil2, 4 * eil, 0, -6 * eil2, 2 * eil],
                       [-eal, 0, 0, eal, 0, 0],
                       [0, -12 * eil3, -6 * eil2, 0, 12 * eil3, -6 * eil2],
                       [0, 6 * eil2, 2 * eil, 0, -6 * eil2, 4 * eil]])
    r = np.array([[cx, cy, 0, 0, 0, 0],
                  [-cy, cx, 0, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, cx, cy, 0],
                  [0, 0, 0, -cy, cx, 0],
                  [0, 0, 0, 0, 0, 1]])
    kglobal = np.transpose(r) @ klocal @ r
    return kglobal


# Trouvage du nb de noeuds
def noeud_max(liste):
    gros_noeud = 0
    for i in range(len(liste)):
        if liste[i][0] > gros_noeud:
            gros_noeud = liste[i][0]
        if liste[i][1] > gros_noeud:
            gros_noeud = liste[i][1]
    return gros_noeud
