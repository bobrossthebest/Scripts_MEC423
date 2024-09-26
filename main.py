# imports externes

import numpy as np

# imports internes

import Modules.Treillis_2D as T2d, Modules.Treillis_3D as T3d, Modules.Poutre_1D as P1d, Modules.Charpente_2D as C2d
import Modules.Surfacique as Sq, Modules.Jeux as Jx, Modules.Ressorts as Rt, Modules.Appui_incliné as App
import Modules.Mécanisme as Mcn, Modules.Flambage as Flb, Modules.Thermique as Thm


# Menu est la fonction qui interagit avec l'utilsateur
# L'affichage est aussi géré par menu


def menu():
    np.set_printoptions(precision=4, linewidth=int(np.inf))

    pass


if __name__ == '__main__': menu()
