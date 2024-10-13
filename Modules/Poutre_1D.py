# Unites : N-mm-MPa

import numpy as np
import math

from Modules.Fonctions_partagées import (calculer_k_poutre1d, assembler_matrice, extraire_matrice, extraire_vecteur,
                                         reconstruire_vecteur, assembler_vecteur)
from Modules.Treillis_3D import elements

# Unites :
F = input("\nQuelle est l'unité de mesure de force?\t\t")
L = input("Quelle est l'unité de mesure de longueur?\t")
M = input("Quelle est l'unité de mesure du moment?\t")
P = input("Quelle est l'unité de mesure de contrainte?\t")
ddl_par_ele = 4

# ---------
# Fonctions
# ---------

# Calcul du vecteur des charges equivalentes dues a une charge repartie constante appliquee sur un element Poutre1D
def calculer_feq_poutre1d(q, l_poutre):
    feq = np.array([[q * l_poutre / 2],
                    [q * l_poutre ** 2 / 12],
                    [q * l_poutre / 2],
                    [-q * l_poutre ** 2 / 12]])
    return feq


# Calcul par interpolation du deplacement transversal dans un element Poutre1D
def calculer_deplacement_poutre1d(u_tot, ddl, l_poutre, x):
    vi = u_tot[ddl[0] - 1][0]
    ti = u_tot[ddl[1] - 1][0]
    vj = u_tot[ddl[2] - 1][0]
    tj = u_tot[ddl[3] - 1][0]
    dep = ((2 * x ** 3 / l_poutre ** 3 - 3 * x ** 2 / l_poutre ** 2 + 1) * vi +
           (x ** 3 / l_poutre ** 2 - 2 * x ** 2 / l_poutre + x) * ti +
           (-2 * x ** 3 / l_poutre ** 3 + 3 * x ** 2 / l_poutre ** 2) * vj +
           (x ** 3 / l_poutre ** 2 - x ** 2 / l_poutre) * tj)
    return dep


# Calcul de la contrainte dans un element Poutre1D
def calculer_contrainte_poutre1d(u_tot, ddl, e, l_poutre, x, y):
    vi = u_tot[ddl[0] - 1][0]
    ti = Utot[ddl[1] - 1][0]
    vj = Utot[ddl[2] - 1][0]
    tj = Utot[ddl[3] - 1][0]
    sigma = -y * e * ((12 * x / l_poutre ** 3 - 6 / l_poutre ** 2) * vi +
                      (6 * x / l_poutre ** 2 - 4 / l_poutre) * ti +
                      (-12 * x / l_poutre ** 3 + 6 / l_poutre ** 2) * vj +
                      (6 * x / l_poutre ** 2 - 2 / l_poutre) * tj)
    return sigma

# Calcul de l'inertie en Iz pour les calculs
def calcul_Iz():

    while True:
        type_element = print(f"Quelle type de poutre avez-vouz:\trectangle\ttriangle\tcercle\tdemi-cercle\tcercle-mince\n"
                           "Écrivez exactement la réponse comme elle est écrite.")
        if type_element == "rectangle":
            b = float(print(f"Quelle est la valeur de la base?"))
            h = float(print(f"Quelle est la valeur de la hauteur?"))
            Iz = (b*(h**3))/12
            return Iz
        elif type_element == "triangle":
            b = float(print(f"Quelle est la valeur de la base?"))
            h = float(print(f"Quelle est la valeur de la hauteur?"))
            Iz = (b*(h**3))/36
            return Iz
        elif type_element == "cercle":
            r = float(print(f"Quelle est la valeur du rayon?"))
            Iz = (math.pi*(r**4))/4
            return Iz
        elif type_element == "demi-cercle":
            r = float(print(f"Quelle est la valeur du rayon?"))
            Iz = math.pi*(r**4)*((1/8)-(8/9*math.pi**2))
            return Iz
        elif type_element == "cercle-mince":
            rm = float(print(f"Quelle est la valeur du rayon moyen?"))
            t = float(print(f"Quelle est la valeur de l'epaisseur?"))
            Iz = math.pi*(rm**3)*t
            return Iz
        else:
            print('Erreur de saisie, veuillez recommencer.')
            continue



# ----------------------------
# Proprietes de chaque element
# ----------------------------

nb_element = int(input("\nCombien de d'element contient la structure?\t"))

#Creation d'un dictionnaire avez toutes les cases pour chaque noeuds
elements = {'ddl': [0,0,0,0] * nb_element, 'L': [0] * nb_element, 'Iz': [0] * nb_element,
              'E': [0] * nb_element, 'ymax': [0] * nb_element, 'k': [0] * nb_element, 'feq': [0]*nb_element }

for element in range(nb_element):
    while True:
        try:
            for i in range(len(elements['ddl'][element])):
                elements['ddl'][element][i] = int(input(f"l'élement {element+1} quel est le {i+1} degrée de liberté?"))

        except(ValueError, SyntaxError, TypeError):
            continue
        break
    while True:
        try:
           elements['L'][element] = float(input(f"Quelle est la longueur de l'élément {element+1}?"))
        except(ValueError, SyntaxError, TypeError):
            continue
        break
    while True:
        try:
           elements['L'][element] = float(input(f"Quelle est la longueur de l'élément {element+1}?"))
        except(ValueError, SyntaxError, TypeError):
            continue
        break

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
