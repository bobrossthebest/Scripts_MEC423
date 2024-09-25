# Unites : lbf-po-psi

import numpy as np
np.set_printoptions(precision=4,linewidth=np.inf)

# ---------
# Fonctions
# ---------

# Calcul de la matrice de rigidite d'un element Barre3D
def calculer_k_barre3d(E,A,xi,yi,zi,xj,yj,zj):
    L = ((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)**0.5
    cx = (xj-xi)/L
    cy = (yj-yi)/L
    cz = (zj-zi)/L
    k = (E*A/L)*np.array([[ cx**2,  cx*cy,  cx*cz, -cx**2, -cx*cy, -cx*cz],
                          [ cx*cy,  cy**2,  cy*cz, -cx*cy, -cy**2, -cy*cz],
                          [ cx*cz,  cy*cz,  cz**2, -cx*cz, -cy*cz, -cz**2],
                          [-cx**2, -cx*cy, -cx*cz,  cx**2,  cx*cy,  cx*cz],
                          [-cx*cy, -cy**2, -cy*cz,  cx*cy,  cy**2,  cy*cz],
                          [-cx*cz, -cy*cz, -cz**2,  cx*cz,  cy*cz,  cz**2]])
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

# Calcul de la contrainte dans un element Barre3D
def calculer_contrainte_barre3d(Utot,ddl,E,alpha,dT,xi,yi,zi,xj,yj,zj):
    Ui = Utot[ddl[0]-1][0]
    Vi = Utot[ddl[1]-1][0]
    Wi = Utot[ddl[2]-1][0]
    Uj = Utot[ddl[3]-1][0]
    Vj = Utot[ddl[4]-1][0]
    Wj = Utot[ddl[5]-1][0]
    L = ((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)**0.5
    cx = (xj-xi)/L
    cy = (yj-yi)/L
    cz = (zj-zi)/L
    sigma = (E/L)*np.array([[-cx,-cy,-cz,cx,cy,cz]])@np.array([[Ui],[Vi],[Wi],[Uj],[Vj],[Wj]]) - np.array([[E*alpha*dT]])
    return sigma[0][0]

# ----------------------------
# Proprietes de chaque element
# ----------------------------

ddl1          = np.array([1,2,3,10,11,12])
xi1, yi1, zi1 =  0, 60, 84
xj1, yj1, zj1 = 24, 60, 48
E1            = 30e6
A1            = 0.2
k1  = calculer_k_barre3d(E1,A1,xi1,yi1,zi1,xj1,yj1,zj1)

ddl2          = np.array([4,5,6,10,11,12])
xi2, yi2, zi2 =  0, 84, 48
xj2, yj2, zj2 = 24, 60, 48
E2            = 30e6
A2            = 0.2
k2  = calculer_k_barre3d(E2,A2,xi2,yi2,zi2,xj2,yj2,zj2)

ddl3          = np.array([7,8,9,10,11,12])
xi3, yi3, zi3 = 24, 60,  0
xj3, yj3, zj3 = 24, 60, 48
E3            = 30e6
A3            = 0.2
k3  = calculer_k_barre3d(E3,A3,xi3,yi3,zi3,xj3,yj3,zj3)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((12,12))
Ktot = assembler_matrice(Ktot,k1,ddl1,ddl1)
Ktot = assembler_matrice(Ktot,k2,ddl2,ddl2)
Ktot = assembler_matrice(Ktot,k3,ddl3,ddl3)

# -------------------------
# Conditions aux frontieres
# -------------------------

ddlUc = np.array([ 1,  2,  3,  4,  5,  6,  7,  8,  9,    11])
Uc    = np.array([[0],[0],[0],[0],[0],[0],[0],[0],[0],[-1.2]])
ddlFc = np.array([10, 12])
Fc    = np.array([[0],[0]])

# ---------------
# Partitionnement
# ---------------

Kic = extraire_matrice(Ktot,ddlFc,ddlFc)
Kcc = extraire_matrice(Ktot,ddlFc,ddlUc)
Kii = extraire_matrice(Ktot,ddlUc,ddlFc)
Kci = extraire_matrice(Ktot,ddlUc,ddlUc)

# --------
# Solution
# --------

Ui = np.linalg.inv(Kic)@(Fc - Kcc@Uc)
Fi = Kii@Ui + Kci@Uc

# --------------
# Reconstruction
# --------------

Utot = reconstruire_vecteur(Uc,ddlUc,Ui,ddlFc)
Ftot = reconstruire_vecteur(Fc,ddlFc,Fi,ddlUc)

# --------
# Reponses
# --------

Fy_F = extraire_vecteur(Ftot,[11])
print('Fy_F = %.0f lbf' % Fy_F[0][0])

sigma3 = calculer_contrainte_barre3d(Utot,ddl3,E3,0,0,xi3,yi3,zi3,xj3,yj3,zj3)
print('sigma_FD = %.0f psi' % sigma3)