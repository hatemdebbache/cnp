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
gALPHA = 2.
gOMEGA2 = 1.
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo3()
# ==============================================================
def exo1():
    a = -1.
    b = 1.
    ya = y_exacte(a)
    c, t, yE = edo.Euler1(Ali, a, b, ya, 20, 40)
    c, t, yH = edo.Heun1(Ali, a, b, ya, 20, 2)
    c, t, yP = edo.PointMilieu1(Ali, a, b, ya, 20, 2)
    c, t, yR = edo.RungeKutta1(Ali, a, b, ya, 20, 1)
    
    graphe.FixeEchelle(-1., 1., 0., 1.1)
    graphe.TraceAxes()
    graphe.TraceFonc(y_exacte, -1., 1., epaisseur=3.)
    graphe.TracePoints(t, yE, couleur='blue', epaisseur=5.)
    # graphe.TracePoints(t, yH, couleur='white', epaisseur=4.)
    # graphe.TracePoints(t, yP, couleur='cyan', epaisseur=4.)
    graphe.TracePoints(t, yR, couleur='red', epaisseur=5.)
    plt.show()
# ==============================================================
def exo2():
    c, t, x, v = edo.RungeKutta2(ressort, 0., 4*pi, 1., 0., 100, 1000)
    graphe.FixeEchelle(-1.1, 1.1, -1.1, 1.1)
    # graphe.FixeEchelle(0., 4*pi, -1., 1.)
    graphe.TraceAxes()
    # graphe.TracePoints(t, x, couleur='green', relie=True, epaisseur=3.)
    # graphe.TracePoints(t, v, couleur='pink', relie=True, epaisseur=3.)
    # graphe.TraceFonc(exacte, 0., pi, couleur = 'red', epaisseur=2.)
    graphe.TracePoints(x, v, couleur='cyan', epaisseur=5.)
    plt.show()
# ==============================================================
def exo3():
    a = -1.
    b = 1.
    ya = y_exacte(a)
    c, t, y = edo.RungeKutta1(Ali, a, b, ya, 20, 1)
    vr = y[10]
    er = abs(vr - 1.)
    
    c, t, y = edo.PointMilieu1(Ali, a, b, ya, 20, 1)
    v1 = y[10]
    e1 = abs(v1 - 1.)
    
    c, t, y = edo.PointMilieu1(Ali, a, b, ya, 20, 2)
    v2 = y[10]
    e2 = abs(v2 - 1.)
    
    v3 = v2 + (v2 - v1)/3.
    e3 = abs(v3 - 1.)
    
    print(e1)
    print(e2)
    print(e3)
    print(er)
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
    v_point = -gALPHA*v - gOMEGA2*x
    return x_point, v_point
# ==============================================================
def exacte(t):
    return cos(sqrt(gOMEGA2)*t)
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