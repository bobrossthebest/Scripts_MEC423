# Unités : N-mm-MPa

import numpy as np
np.set_printoptions(precision=4,linewidth=np.inf)

# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Barre2D
def calculer_k_barre2d(E,A,xi,yi,xj,yj):
    L = ((xj-xi)**2+(yj-yi)**2)**0.5
    cx = (xj-xi)/L
    cy = (yj-yi)/L
    k = (E*A/L)*np.array([[ cx**2,  cx*cy, -cx**2, -cx*cy],
                          [ cx*cy,  cy**2, -cx*cy, -cy**2],
                          [-cx**2, -cx*cy,  cx**2,  cx*cy],
                          [-cx*cy, -cy**2,  cx*cy,  cy**2]])
    return k

# Assemblage d'une matrice
def assembler_matrice(Ktot,k,lig,col):
    for i in range(len(lig)):
        for j in range(len(col)):
            Ktot[lig[i]-1][col[j]-1] += k[i][j]
    return Ktot

# Extraction d'une sous-matrice
def extraire_matrice(Kin,lig,col):
    Kout = np.zeros((len(lig),len(col)))
    for i in range(len(lig)):
        for j in range(len(col)):
            Kout[i][j] = Kin[lig[i]-1][col[j]-1]
    return Kout

# Reconstruction d'un vecteur a partir de deux vecteurs
def reconstruire_vecteur(Vin1,lig1,Vin2,lig2):
    Vout = np.zeros((len(lig1)+len(lig2),1))
    for i in range(len(lig1)):
        Vout[lig1[i]-1][0] = Vin1[i][0]
    for i in range(len(lig2)):
        Vout[lig2[i]-1][0] = Vin2[i][0]
    return Vout

# Extraction d'un vecteur
def extraire_vecteur(Vin,lig):
    Vout = np.zeros((len(lig),1))
    for i in range(len(lig)):
        Vout[i][0] = Vin[lig[i]-1][0]
    return Vout

# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1     = np.array([1,2,3,4])
xi1, yi1 =   0,   0
xj1, yj1 = 800, 600
E1       = 2e5
A1       = 24
k1 = calculer_k_barre2d(E1,A1,xi1,yi1,xj1,yj1)

ddl2     = np.array([3,4,5,6])
xi2, yi2 = 800, 600
xj2, yj2 = 800,   0
E2       = 2e5
A2       = 24
k2 = calculer_k_barre2d(E2,A2,xi2,yi2,xj2,yj2)

ddl3     = np.array([3,4,7,8])
xi3, yi3 =  800, 600
xj3, yj3 = 1600,   0
E3       = 2e5
A3       = 24
k3 = calculer_k_barre2d(E3,A3,xi3,yi3,xj3,yj3)

ddl4     = np.array([1,2,5,6])
xi4, yi4 =   0,   0
xj4, yj4 = 800,   0
E4       = 2e5
A4       = 24
k4 = calculer_k_barre2d(E4,A4,xi4,yi4,xj4,yj4)

ddl5     = np.array([5,6,7,8])
xi5, yi5 =  800,   0
xj5, yj5 = 1600,   0
E5       = 2e5
A5       = 24
k5 = calculer_k_barre2d(E5,A5,xi5,yi5,xj5,yj5)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((8,8))
Ktot = assembler_matrice(Ktot,k1,ddl1,ddl1)
Ktot = assembler_matrice(Ktot,k2,ddl2,ddl2)
Ktot = assembler_matrice(Ktot,k3,ddl3,ddl3)
Ktot = assembler_matrice(Ktot,k4,ddl4,ddl4)
Ktot = assembler_matrice(Ktot,k5,ddl5,ddl5)

# --------
# Couplage
# --------

matA = np.array([[0,0,0,0,0,0,1,1]])

# ------------
# Augmentation
# ------------

Kaug = np.zeros((9,9))
ddlinit = np.array([1,2,3,4,5,6,7,8])
ddlaugm = np.array([9])
Kaug = assembler_matrice(Kaug,Ktot,ddlinit,ddlinit)
Kaug = assembler_matrice(Kaug,matA,ddlaugm,ddlinit)
Kaug = assembler_matrice(Kaug,matA.transpose(),ddlinit,ddlaugm)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([ 1,  2,  6])
Uc    = np.array([[0],[0],[0]])
ddlFc = np.array([    3,  4,  5,  7,  8,  9])
Fc    = np.array([[1000],[0],[0],[0],[0],[0]])

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Kaug,ddlFc,ddlFc)
Kcc = extraire_matrice(Kaug,ddlFc,ddlUc)
Kii = extraire_matrice(Kaug,ddlUc,ddlFc)
Kci = extraire_matrice(Kaug,ddlUc,ddlUc)

# --------
# Solution
# --------

Ui = np.linalg.inv(Kic)@(Fc - Kcc@Uc)
Fi = Kii@Ui + Kci@Uc

# --------------
# Reconstruction
# --------------

# Incluant les multiplicateurs de Lagrange
Utot = reconstruire_vecteur(Uc,ddlUc,Ui,ddlFc)

# Retrait des multiplicateurs de Lagrange
Utot = extraire_vecteur(Utot,ddlinit)

# Obtention des forces externes
Ftot = Ktot @ Utot

# --------
# Reponses
# --------

Ux_F = extraire_vecteur(Utot,[7])
Uy_F = extraire_vecteur(Utot,[8])
Rx_F = extraire_vecteur(Ftot,[7])
Ry_F = extraire_vecteur(Ftot,[8])
print('Deplacement horizontal = %.3f mm' % Ux_F[0][0])
print('Deplacement vectical   = %.3f mm' % Uy_F[0][0])
print('Reaction horizontale   = %.1f N' % Rx_F[0][0])
print('Reaction vecticale     = %.1f N' % Ry_F[0][0])