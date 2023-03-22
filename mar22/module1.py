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
    exo2()
# =============================================================
def exo1():
    graphe.FixeEchelle(-10., 10., -10., 10.)
    graphe.TraceAxes()
    c, t, y, z = edo.RungeKutta2(samia, 0., 2., -0.5, 0.3, 1000, 100)
    graphe.TracePoints(y, z)
    c, t, y, z = edo.RungeKutta2(samia, 2., 0., -0.5, 0.3, 1000, 100)
    graphe.TracePoints(y, z)
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 10., 0.2, 0.2, 1000, 100)
    graphe.TracePoints(y, z, couleur='blue')
    c, t, y, z = edo.RungeKutta2(samia, 10., 0., 0.2, 0.2, 1000, 100)
    graphe.TracePoints(y, z, couleur='blue')
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 2., -0.1, -0.1, 1000, 100)
    graphe.TracePoints(y, z, couleur='green')
    c, t, y, z = edo.RungeKutta2(samia, 2., 0., -0.1, -0.1, 1000, 100)
    graphe.TracePoints(y, z, couleur='green')
    
    plt.show()
# ==============================================================
def exo2():
    graphe.FixeEchelle(-10., 10., -10., 10.)
    graphe.TraceAxes()
    c, t, y, z = edo.RungeKutta2(samia2, 0., 50., -0.5, 0.3, 1000, 100)
    graphe.TracePoints(y, z)
    c, t, y, z = edo.RungeKutta2(samia2, 50., 0., -0.5, 0.3, 1000, 100)
    graphe.TracePoints(y, z)
    
    c, t, y, z = edo.RungeKutta2(samia2, 0., 50., 0.5, 0.3, 1000, 100)
    graphe.TracePoints(y, z, couleur='blue')
    c, t, y, z = edo.RungeKutta2(samia2, 50., 0., 0.5, 0.3, 1000, 100)
    graphe.TracePoints(y, z, couleur='blue')
    
    c, t, y, z = edo.RungeKutta2(samia2, 0., 50., -1., -2, 1000, 100)
    graphe.TracePoints(y, z, couleur='green')
    c, t, y, z = edo.RungeKutta2(samia2, 50., 0., -1., -2, 1000, 100)
    graphe.TracePoints(y, z, couleur='green')
    
    plt.show()
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
    return 3.*y + 4.*z, 5.*y + 2.*z
# ==============================================================
def samia2(t, y, z):
    return -y + z, (-2. + gALPHA)*y - z
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================