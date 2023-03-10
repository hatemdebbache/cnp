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
import graphe
import edo
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
gALPHA = 10.
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
    a = -1.
    b = 1.
    ya = y_exacte(a)
    c, t, y = edo.Euler1(Ali, a, b, ya, 20, 2000)
    print(c, abs(1.-y[10]))
    graphe.FixeEchelle(-1., 1., 0., 1.1)
    graphe.TraceAxes()
    graphe.TraceFonc(y_exacte, -1., 1.)
    graphe.TracePoints(t, y)
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
def exo5():
    pass
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def Ali(t, y):
    return -2.*gALPHA*t*y*y
# ==============================================================
def y_exacte(t):
    return 1./(1. + gALPHA*t*t)
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================