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
import cbase
import edo
import graphe

# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
cs = 4         # Numbre of most significant digits to be printed
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    # exo1()
    exo2()
    # exo3()
    # exo4()
    # exo5()
# ==============================================================
def exo1():
    print("==================== EXO 1 =====================")
    # # Données
    # print("C = 0\n") ; c = 0.
    # def eq_diff_1(t,y,z):
    #     return c*y + z, -(1. + 10*gBETA)*y + c*z
    # gOMEGA = sqrt(1 + 10*gBETA)
    # pc = [0., 0.]
    # to = 0. ; tf = 1.5
    # ya = 5*gBETA ; za = 0.
    # h = 0.01 ; npts = int((tf - to)/h)

    # # SOL EXACTE:
    # A = 5*gBETA ; B = 0.
    # def sol_exacte_1(t):
    #     return A*cos(gOMEGA*t) + B*sin(gOMEGA*t)
    # V = sol_exacte_1(1.5)
    # print(f'A = {A:.4g} \tB = {B:.4g}\ny(t=1.5s) = {V:.9g}\n')
    
    # Nsub = 150 ; m = npts//Nsub
    # c, t, y, z = edo.RungeKutta2(eq_diff_1, to, tf, ya, za, Nsub, m)
    # Wrk = y[-1] ; Erk = abs(Wrk - V)
    # print(f'Wrk = {Wrk:.9g} \tErk = {Erk:.2g}\n')    
    
    # # Graphe
    # # ymin = -2. ; ymax = 2.
    # # zmin = -2. ; zmax = 2.
    # # graphe.FixeEchelle(ymin, ymax, zmin, zmax)
    # # graphe.TraceAxes()
    # # graphe.TracePoints(y, z, couleur='green', relie=True)
    # # graphe.TracePoints(t, z, relie=True)
    # # plt.show()
    
    # # Qst 6
    # # i_min = np.argmin(z) ; t_etoile = t[i_min]
    # # print(f't_etoile = {t_etoile:.3g}\n')
    
    
    # print("C != 0\n") ; c = -2.
    # to = 0. ; tf = -0.8
    # dt = 0.001 ; npts = int(abs((tf - to)/dt))
    # Nsub = 500 ; m = npts//Nsub

    # ymin = -4. ; ymax = 4.
    # zmin = -4. ; zmax = 4.
    # graphe.FixeEchelle(ymin, ymax, zmin, zmax)
    # graphe.TraceAxes()
    # r = 0.001 ; theta = np.linspace(0,2*pi,num=10,endpoint=False)
    # Yo = pc[0] + r*np.cos(theta)
    # Zo = pc[0] + r*np.sin(theta)
    # for ya, za in zip(Yo, Zo):
    #     c, t, y, z = edo.RungeKutta2(eq_diff_1, to, tf, ya, za, Nsub, m)
    #     graphe.TracePoints(y, z, relie=True)
    
    # # Qst 2
    # c, t, y7, z7 = edo.RungeKutta2(eq_diff_1, to, tf, Yo[7], Zo[7], Nsub, m)
    # cs = 4
    # print(f'\ny7(100) = {y7[100]:.{cs}g} \tz7(100) = {z7[100]:.{cs}g}\n')
    
    # # Qst 3
    # c, t, y5, z5 = edo.RungeKutta2(eq_diff_1, to, tf, Yo[5], Zo[5], Nsub, m)
    # cs = 4
    # print(f'\ny5(100) = {y5[100]:.{cs}g} \tz5(100) = {z5[100]:.{cs}g}\n')
    
    # Qst 4 
    c = -2.1
    A = np.array([[c, 1.], [-(1 + 10*gBETA), c]])
    xmin = -5. ; xmax = 4.
    ymin = -2. ; ymax = 10.
    graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    graphe.TraceAxes()
    sigma = np.linspace(xmin, xmax) ; fct = (0.5*sigma)**2
    graphe.TracePoints(sigma, fct, couleur='white')
    plt.axhline(c*c + 1 + 10*gBETA, lw=1., color='white')
    plt.axvline(2*c, lw=1., color='white')
    print("================================================\n")
# ==============================================================
def exo2():
    print("==================== EXO 2 =====================")
    
    def eq_diff_2(t, x, y):
        return x*x + c, -10*gBETA*y
        
    print("================================================\n")
# ==============================================================
def exo3():
    print("==================== EXO 3 =====================")
    a = 1
    b = 1
    c = 1
    d = gBETA
    dt = 0.1
    ya = 5. ; za = 1.
    def eq_diff_3(t, y, z):
        return a*y*(1 - z/b) , -c*z*(1 - y/d)
    c, t, y, z = edo.RungeKutta2(eq_diff_3, 0., 100, ya, za, 1000, 1)
    
    print("================================================\n")
# ==============================================================
def exo4():
    print("==================== EXO 4 =====================")
    
    
    print("================================================\n")
# ==============================================================
def exo5():
    print("==================== EXO 5 =====================")
    
    
    print("================================================\n")
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def f(x):
    pass
# ==============================================================
def g(x):
    pass
# ==============================================================

# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================