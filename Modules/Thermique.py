# Unités : m-W-degC

import numpy as np, sympy as sp

from Modules.Fonctions_partagées import (assembler_matrice, assembler_vecteur)


sp.var('x y')


# ---------
# Fonctions
# ---------

# Calcul des fonctions d'interpolation et de leurs dérivées d'un élément RECT4
def calculer_deriver_fonctions_interpolation_rect4(xmin, xmax, ymin, ymax):
    ni = ((xmax - x) * (ymax - y)) / ((xmax - xmin) * (ymax - ymin))
    nj = ((x - xmin) * (ymax - y)) / ((xmax - xmin) * (ymax - ymin))
    nk = ((x - xmin) * (y - ymin)) / ((xmax - xmin) * (ymax - ymin))
    nl = ((xmax - x) * (y - ymin)) / ((xmax - xmin) * (ymax - ymin))
    n = sp.Matrix([ni, nj, nk, nl])
    dn_dx = sp.Matrix([ni.diff(x), nj.diff(x), nk.diff(x), nl.diff(x)])
    dn_dy = sp.Matrix([ni.diff(y), nj.diff(y), nk.diff(y), nl.diff(y)])
    return n, dn_dx, dn_dy


# ----------------------------
# Proprietes de chaque element
# ----------------------------

# Parametres communs aux elements et conditions frontieres thermiques
t = 0.01
kT = 22
h = 50
Tf = 36
q14 = 30 / (0.5 * t)

# Element 1
ddl1 = np.array([1, 2, 5, 4])
xmin1 = 0
xmax1 = 0.25
ymin1 = 0
ymax1 = 0.25
N1, dNdx1, dNdy1 = calculer_deriver_fonctions_interpolation_rect4(xmin1, xmax1, ymin1, ymax1)
Kk1 = t * kT * sp.integrate(sp.integrate(dNdx1 * dNdx1.T + dNdy1 * dNdy1.T, (x, xmin1, xmax1)), (y, ymin1, ymax1))
Kh1 = t * h * (sp.integrate((N1 * N1.T).subs(y, ymin1), (x, xmin1, xmax1)))
Fh1 = t * h * Tf * sp.integrate(N1.subs(y, ymin1), (x, xmin1, xmax1))
Fq1 = t * q14 * sp.integrate(N1.subs(x, xmin1), (y, ymin1, ymax1))
K1 = np.array(Kk1 + Kh1)
F1 = np.array(Fh1 + Fq1)

# Element 2
ddl2 = np.array([2, 3, 6, 5])
xmin2 = 0.25
xmax2 = 0.5
ymin2 = 0
ymax2 = 0.25
N2, dNdx2, dNdy2 = calculer_deriver_fonctions_interpolation_rect4(xmin2, xmax2, ymin2, ymax2)
Kk2 = t * kT * sp.integrate(sp.integrate(dNdx2 * dNdx2.T + dNdy2 * dNdy2.T, (x, xmin2, xmax2)), (y, ymin2, ymax2))
Kh2 = t * h * (sp.integrate((N2 * N2.T).subs(y, ymin2), (x, xmin2, xmax2))) + t * h * (
    sp.integrate((N2 * N2.T).subs(x, xmax2), (y, ymin2, ymax2)))
Fh2 = (t * h * Tf * sp.integrate(N2.subs(y, ymin2), (x, xmin2, xmax2)) +
       t * h * Tf * sp.integrate(N2.subs(x, xmax2), (y, ymin2, ymax2)))
K2 = np.array(Kk2 + Kh2)
F2 = np.array(Fh2)

# ----------
# Assemblage
# ----------

Ktot = np.zeros((6, 6))
Ktot = assembler_matrice(Ktot, K1, ddl1, ddl1)
Ktot = assembler_matrice(Ktot, K2, ddl2, ddl2)

Ftot = np.zeros((6, 1))
Ftot = assembler_vecteur(Ftot, F1, ddl1)
Ftot = assembler_vecteur(Ftot, F2, ddl2)

# --------
# Solution
# --------

Ttot = np.linalg.inv(Ktot) @ Ftot

# --------
# Reponses
# --------

print('Tmax = %.1f C' % max(Ttot)[0])
print('Tmin =  %.1f C' % min(Ttot)[0])
