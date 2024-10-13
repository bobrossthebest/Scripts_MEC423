# Unites : N-mm-MPa

import numpy as np
import math
import keyboard

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
        try:
            type_element = print(f"Quelle type de poutre avez-vouz:\trectangle\ttriangle\tcercle\tdemi-cercle\tcercle-mince\n"
                           "Écrivez exactement la réponse comme elle est écrite.")
        except(ValueError, SyntaxError, TypeError):
            continue

        if type_element == "rectangle":
            b = float(print(f"Quelle est la valeur de la base en {L}?"))
            h = float(print(f"Quelle est la valeur de la hauteur en {L}?"))
            Iz = (b*(h**3))/12
            return Iz
        elif type_element == "triangle":
            b = float(print(f"Quelle est la valeur de la base en {L}?"))
            h = float(print(f"Quelle est la valeur de la hauteur en {L}?"))
            Iz = (b*(h**3))/36
            return Iz
        elif type_element == "cercle":
            r = float(print(f"Quelle est la valeur du rayon en {L}?"))
            Iz = (math.pi*(r**4))/4
            return Iz
        elif type_element == "demi-cercle":
            r = float(print(f"Quelle est la valeur du rayon en {L}?"))
            Iz = math.pi*(r**4)*((1/8)-(8/9*math.pi**2))
            return Iz
        elif type_element == "cercle-mince":
            rm = float(print(f"Quelle est la valeur du rayon moyen en {L}?"))
            t = float(print(f"Quelle est la valeur de l'epaisseur en {L}?"))
            Iz = math.pi*(rm**3)*t
            return Iz
        else:
            print('Erreur de saisie, veuillez recommencer.')
            continue



# ----------------------------
# Proprietes de chaque element
# ----------------------------

nb_element = int(input("\nCombien d'element contient la structure?\t"))

#Creation d'un dictionnaire avez toutes les cases pour chaque noeuds
elements = {'ddl': [0,0,0,0] * nb_element, 'L': [0] * nb_element, 'Iz': [0] * nb_element,
              'E': [0] * nb_element,'q': [0] * nb_element, 'ymax': [0] * nb_element, 'k': [0] * nb_element, 'feq': [0]*nb_element }

redo = True
while redo is True:
    for element in range(nb_element):
        #Te demande les ddl de la membrure
        while True:
            try:
                for i in range(len(elements['ddl'][element])):
                    elements['ddl'][element][i] = int(input(f"l'élement {element+1} quel est le {i+1} degrée de liberté?"))

            except(ValueError, SyntaxError, TypeError):
                continue
            break
        #Te demande la longueur de l'element
        while True:
            try:
               elements['L'][element] = float(input(f"Quelle est la longueur de l'élément {element+1} en {L}?"))
            except(ValueError, SyntaxError, TypeError):
                continue
            break
        #T'envoie dans la function Iz pour calculer l'inertie
        elements['Iz'][element] = calcul_Iz()
        #Demande le module d'elasticit/ de l'element
        while True:
            try:
               elements['E'][element] = float(input(f"Quelle est le module d'elasticite en {P} de l'élément {element}?"))
            except(ValueError, SyntaxError, TypeError):
                continue
            break
        #Demande si il y a une charge repartie
        print(f"Appuyez sur '1' si l'élément {element+1} a une charge répartie, ou '0' pour pas de charge.")
        while True:
            if keyboard.is_pressed('1'):
                q_present = 1
                print("Vous avez appuyé sur 1.")
                while True:
                    try:
                        elements['q'][element] = float(input(f"Quelle est la charge répartie en {F}/{L} de l'élément {element}?"))
                    except(ValueError, SyntaxError, TypeError):
                        continue
                    break
                break
            elif keyboard.is_pressed('0'):
                q_present = 0
                print("Vous avez appuyé sur 0.")
                elements['q'][element] = False
                break
        #Demande la fibre la plus éloigné de la fibre neutre
        while True:
            try:
               elements['ymax'][element] = float(input(f"Quelle est le ymax en {L} de l'élément {element}?"))
            except(ValueError, SyntaxError, TypeError):
                continue
            break
        #Calcul le k de l'element
        elements['k'][element] = calculer_k_poutre1d(elements['E'][element],elements['Iz'][element],elements['L'][element])
        #Si charge repartie est presente calcul de feq
        if q_present == 1:
            elements['feq'][element] = calculer_feq_poutre1d(elements['q'][element],elements['L'][element])
        elif q_present == 0:
            elements['feq'][element] = False
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))


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
# Le plus 1 donne le nombre de noeud et *2 car 2 degree de liberter par noeud
Ktot = np.zeros(((nb_element+1)*2,((nb_element+1)*2 )))
for i in range(nb_element):
    Ktot = assembler_matrice(Ktot, elements['k'][i], elements['ddl'][i], elements['ddl'][i])

Feqtot = np.zeros((((nb_element+1)*2, 1)))
for i in range(nb_element):
    if elements['feq'][i] == False:
        continue
    else:
        Feqtot = assembler_vecteur(Feqtot, elements['feq'][i], elements['ddl'][i])
# -------------------------
# Conditions aux frontieres
# -------------------------
redo = True
while redo is True:
    while True:
        try:
            nb_Uc = int(input('Combien de déplacements sont connus? '))
        except (ValueError, SyntaxError, TypeError):
            continue
        break
    ddlUc = [0] * nb_Uc
    Uc = np.zeros((nb_Uc, 1))
    for i in range(nb_Uc):
        while True:
            try:
                ddlUc[i] = int(input(f'Numéro du ddl déplacement connu #{i + 1}: '))
                Uc[i][0] = eval(input(f'Déplacement en {L} de U{ddlUc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    for i in range(nb_Uc):
        print(f'U{ddlUc[i]}:\t{Uc[i][0]} {L}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

redo = True
while redo is True:
    try:
        nb_Fc = int(input('Combien de forces sont connues? '))
    except (ValueError, SyntaxError, TypeError):
        continue
    ddlFc = [0] * nb_Fc
    Fc = np.zeros((nb_Fc, 1))
    for i in range(nb_Fc):
        while True:
            try:
                ddlFc[i] = int(input(f'Numéro du ddl de la force connue #{i + 1}: '))
                Fc[i][0] = eval(input(f'Grandeur en {F} de  F{ddlFc[i]}: '))
            except (ValueError, SyntaxError, TypeError):
                print('Valeur erronnée')
                continue
            break
    for i in range(nb_Fc):
        print(f'F{ddlFc[i]}:\t{Fc[i][0]:.2} {F}')
    redo = bool(input('\nAppuyez sur Enter pour passer à la prochaine étape, entrez 1 pour recommencer\n'))

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

print('Ktot:')
for i in range(len(Ktot)):
    print(Ktot[i])
print('Kic:')
for i in range(len(Kic)):
    print(Kic[i])
print('\nKcc:')
for i in range(len(Kcc)):
    print(Kcc[i])
print('\nKii:')
for i in range(len(Kii)):
    print(Kii[i])
print('\nKci:')
for i in range(len(Kci)):
    print(Kci[i])
print('\n')

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
