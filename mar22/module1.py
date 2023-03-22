# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020
Updated on Sat Feb 02 11:37:45 2023

@author: user
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi, cos, sin, log, sqrt
# ----- User-Defined-Modules Importing
sys.path.insert(0, "./../BIBLIO/")
import graphe
import edo
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
gALPHA = 0.
gOMEGA = 1.
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo1()
# =============================================================
def exo1():
    c, t, x, v = edo.RungeKutta2(ressort, 0., 4*pi, 1., 0., 600, 100)
    graphe.FixeEchelle(-1.1, 1.1, -1.1, 1.1)
    # graphe.FixeEchelle(0., 4*pi, -1., 1.)
    graphe.TraceAxes()
    # graphe.TracePoints(t, x, couleur='green', relie=True, epaisseur=3.)
    # graphe.TracePoints(t, v, couleur='pink', relie=True, epaisseur=3.)
    # graphe.TraceFonc(exacte, 0., 4*pi, couleur = 'red', epaisseur=2., npts=1000)
    graphe.TracePoints(x, v, couleur='cyan', epaisseur=5.)
    c, t, x, v = edo.RungeKutta2(ressort, 0., pi, 0.5, 0.7, 600, 100)
    graphe.TracePoints(x, v, couleur='red', epaisseur=5.)
    plt.show()
# ==============================================================
def exo3():
    pass
# ==============================================================
def exo4():
    pass
# ==============================================================
def exo5():
    pass
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def ressort(t, x, v):
    x_point = v
    v_point = -gALPHA*v - gOMEGA*gOMEGA*x
    return x_point, v_point
# ==============================================================
def exacte(t):
    return cos(gOMEGA*t)
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================