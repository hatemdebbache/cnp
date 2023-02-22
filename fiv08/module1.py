# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020
Updated on Sat Feb 02 11:37:45 2023

@author: user
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi, cos, sin, log
# ----- User-Defined-Modules Importing
sys.path.insert(0, "./../BIBLIO/")
import cbase
import graphe
import edo

# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
cs = 6         # Numbre of most significant digits to be printed
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo1()
    # exo2()
    # exo3()
    # exo4()
    # exo5()
# ==============================================================
def exo1():
    Va = f(0.)
    Vb = f(3.)
    x, y = edo.DFC2_GT(samia, 0., 3., Va, Vb, 60, True, True)
    graphe.FixeEchelle(-0.1, 3.1, -0.1, 0.5)
    graphe.TraceAxes()
    graphe.TraceFonc(f, 0., 3., epaisseur=3.)
    graphe.TracePoints(x, y, epaisseur=5.)
    plt.show()
    
    v = f(1.6)
    w = y[32]
    e = abs(v - w)
    print(e, w)
# ==============================================================
def exo2():
    a, b = 0., 6.
    Va = 0.
    Vb = h(b)
    x, y = edo.DFC2_GT(idris, a, b, Va, Vb, 600, False, True)
    graphe.FixeEchelle(a-0.1, b+0.1, -1.5, 1.5)
    graphe.TraceAxes()
    graphe.TraceFonc(h, a, b, epaisseur=1., npts=1001)
    graphe.TracePoints(x, y, epaisseur=4.)
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
def samia(x):
    P = -2.
    Q = -1.
    R = 0.
    return P, Q, R
# ==============================================================
def idris(x):
    P = 2/(2*x + 1)
    Q = -(2*x+1)**2
    R = 0.
    return P, Q, R
# ==============================================================
def f(x):
    return x*exp(-x)
# ==============================================================
def h(x): # Une solution de y'' = 2y'/(2x+1) - (2x+1)²y
          # y(0) = 0 ; y(2*pi) = -0.2070018496552964
    return cos(x*x + x)
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================