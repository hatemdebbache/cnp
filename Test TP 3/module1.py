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
import syslin
import derive
import quadra
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
cs = 5         # Numbre of most significant digits to be printed
cse = 2        # Significant digits numbers for error values
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    # exo1()
    # exo2()
    exo3()
    # exo4()
    # exo5()

# ==============================================================
def exo1():
    print("==================== EXO 1 =====================")
    # P1
    # Q1
    a, b = 0., 2.
    ya, yb = 0., 9*gBETA/16
    Nsub = 50
    x, y = edo.DFC2_GT(eqdiff1, a, b, ya, yb, Nsub, True, True)
    x1 = int(Nsub*(0.8-a)/(b - a))
    x2 = int(Nsub*(1.2-a)/(b - a))
    print(f"y({x[x1]}) = {y[x1]:.{cs}g} \t y({x[x2]}) = {y[x2]:.{cs}g}\n")
    # Q2
    h = 0.02
    Nsub = int((b - a)/h)
    xx, yy = edo.DFC2_GT(eqdiff1, a, b, ya, yb, Nsub, True, True)
    e = abs(y[x2] - yy[2*x2])
    print(f"y({xx[2*x2]}) = {yy[2*x2]:.{cs}g} \t E = {e:.2g}\n")
    
    # P2
    yb = 37*gBETA/20
    Nsub = 50
    x, y = edo.DFC2_GT(eqdiff1, a, b, ya, yb, Nsub, True, False)
    Nsub = 100
    xx, yy = edo.DFC2_GT(eqdiff1, a, b, ya, yb, Nsub, True, False)
    e = abs(y[-1] - yy[-1])
    print(f"y({x[-1]}) = {yy[-1]:.{cs}g} \t E = {e:.{cse}g}\n") 
    print("================================================\n")
# ==============================================================
def exo2():
    print("==================== EXO 2 =====================")
    # Q1
    print("A = 9.0 \t B = 1.0\n")    
    # Q2    
    a, b = 0., 1.
    ya, yb = 0., 9*gBETA*(exp(gBETA) - exp(-9*gBETA))
    h = 0.01
    Nsub = int((b - a)/h)
    x, y = edo.DFC2_GT(eqdiff2, a, b, ya, yb, Nsub, False, False)
    y_exacte = np.array([sol_exacte2(pt) for pt in x])
    E = abs(y - y_exacte)
    Emin = min(E) ; Emax = max(E) ; Emoy = np.mean(E)
    print(f"Pour DFC2GT:\nEmin = {Emin:.{cse}g}\t Emax = {Emax:.{cse}g}\t Emoy = {Emoy:.{cse}g}\n")
    # Q3
    x, y = edo.DFC4_GT(eqdiff2, a, b, ya, yb, Nsub, False, False)
    E = abs(y - y_exacte)
    Emin = min(E) ; Emax = max(E) ; Emoy = np.mean(E)
    print(f"Pour DFC4GT:\nEmin = {Emin:.{cse}g}\t Emax = {Emax:.{cse}g}\t Emoy = {Emoy:.{cse}g}\n")
    
    # Visualisation des solutions
    # plt.clf()
    # plt.scatter(x, y, color='red')
    # plt.plot(x, y_exacte, color='blue')
    print("================================================\n")
# ==============================================================
def exo3():
    print("==================== EXO 3 =====================")
    n = 50
    tmin, tmax = 0.5, 2.
    npts = 1000
    t = np.linspace(tmin, tmax, npts)
    X = np.empty(npts)
    L = 1. ; D = 4. ; U = 1.
    b = np.zeros(n)
    for i in range(npts):
        b[1] = sin(10*gBETA*t[i]*t[i])
        x = syslin.Gauss_Thomas_LDU(L, D, U, b)
        X[i] = syslin.Norme2(x)/np.sqrt(n)
    print(f"t* = {t[np.argmin(X)]:.{cse}g}\n")
    
    graphe.FixeEchelle(tmin, tmax, 0., 0.1)
    graphe.TraceAxes()
    graphe.TracePoints(t, X)
    plt.show()
    print("================================================\n")
# ==============================================================
def exo4():
    print("==================== EXO 4 =====================")
    M = 1.
    m = 1.
    g = 10.
    alph = 10.
    xo = 0. ; vo = 0.
    Nsub = 100
    tmin = 0. ; tmax = 1. ; dt = (tmax - tmin)/Nsub
    t = np.array([tmin + i*dt for i in range(Nsub + 1)])
    x, v, a = np.zeros((3,Nsub+1))
    x[0] = xo ; v[0] = vo
    i = 0
    while(t[i] < tmax): # Integration avec la methode d'Euler
        a[i] = (-alph*v[i] + m*g)/M
        v[i+1] += v[i] + a[i]*dt
        x[i+1] += x[i] + v[i]*dt
        i += 1    
    t02 = int((0.2 - tmin)*Nsub/(tmax - tmin))
    print(f"x(0.2) = {x[t02]:.{cs}g}\t v(0.2) = {v[t02]:.{cs}g}\t a(0.2) = {a[t02]:.{cs}g}\n")
    
    # Plotting results
    plt.clf()
    plt.plot(t, v, label='vitesse')
    plt.plot(t, x, label='position')
    plt.plot(t, a, label='acceleration')
    plt.legend(fontsize=14.)
    plt.show()
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
def eqdiff1(x):
    p = 1/x
    q = 1/(x*x)
    r = gBETA*x
    return p, q, r
# ==============================================================
def eqdiff2(x):
    p = -8*gBETA
    q = 9*gBETA*gBETA
    r = 0.
    return p, q, r
# ==============================================================
def sol_exacte2(x):
    A = 9. ; B = 1.
    return A*exp(gBETA*x) + B*exp(-9*gBETA*x)
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
   
    main()
    
# ==============================================================