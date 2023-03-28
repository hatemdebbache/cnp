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
gALPHA = 2.9
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
    graphe.FixeEchelle(-2., 2., -2., 2.)
    graphe.TraceAxes()
    for k in range(20):
        q = 2.*pi*k/20.
        ya = 0.01*cos(q)
        za = 0.01*sin(q)
        
        c, t, y, z = edo.RungeKutta2(samia, 0., 2., ya, za, 400, 100)
        graphe.TracePoints(y, z, couleur='red')
    
    
    plt.show()
# ==============================================================
def exo2():
    pass
# ==============================================================
def exo3():
    pass
# ==============================================================
def exo4():
    pass
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def samia(t, y, z):
    return 4.*y + z, -2.*y + z
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================