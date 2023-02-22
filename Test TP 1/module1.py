# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020

@author: user
"""
# %matplotlib auto
# %matplotlib inline
from math import exp, pi, cos, sin, log, sqrt
import numpy as np
import matplotlib.pyplot as plt
import sys
# -----
sys.path.insert(0, "./../BIBLIO/")
import graphe
import disque
import approx
import cbase
import racine

# -----
gBETA = 0.1433
# ==============================================================
def main():
    exo1()
    exo2()
    exo3()
    exo4()
    pass
# ==============================================================
def exo1():
    print("############# EXO 1 ##############")
    x = np.linspace(1,11,11)
    y = np.array([-1.57*(1+gBETA**2), gBETA/10, 1.4, 1.1+(1+gBETA)**2, (2+gBETA)**2, 2*pi-gBETA, 7.85, 9+2*gBETA, 2*(5+2*gBETA), 11*(1+gBETA), 14.])
    a, b = approx.Droite_MC(x, y)
    print(f"a = {a:.4g} \t b = {b:.4g}")
    def f(u):
        return a*u + b
    S = approx.EQM(x, y, f)
    print(f"S = {S:.5g}")

    plt.scatter(x, y, color='red')
    plt.plot(x, a*x+b)
# ==============================================================
def exo2():
    print("############# EXO 2 ##############")
    n, X, Y = disque.LireVec2('Data1.txt')
    print(f"N = {n}")
    
    # for i in range(n):
    #     if X[i] > 10:
    #         print(i)
    #         break
    kmin, kmax = 123, 384
    print(f"k1 = {kmin} \t k2 = {kmax}")
    a, b = approx.Droite_MC(X, Y, kmin, kmax)
    print(f"a = {a:.4g} \t b = {b:.4g}")
    def f(u):
        return a*u + b
    S = approx.EQM(X, Y, f, None, kmin, kmax)
    print(f"S = {S:.5g}")
    
    xmin, xmax = 0, 13
    ymin, ymax = -10, 60
    graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    graphe.TraceAxes()
    graphe.TracePoints(X[123:385], Y[123:385], couleur='blue', epaisseur=2.0)
    graphe.TracePoints(X, Y)
    graphe.TraceFonc(f, xmin, xmax, couleur='black', epaisseur=2.5)
    
    plt.show
# ==============================================================
def exo3():
    print("############# EXO 3 ##############")
    n, X, Y = disque.LireVec2('Data2.txt')
    
    C = approx.Poly_MC(X, Y, 4)
    # C = np.polynomial.Polynomial.fit(X, Y, 4).convert().coef
    for i in range(5):
        print(f"C{i+1} = {C[i]:.4g}")
    
    def f(u):
        return cbase.Horner(C, u)
    S = approx.EQM(X, Y, f)
    print(f"S = {S:.5g}")
    
    def w(x):
        return 10 - x*x
    def d(x):
        return f(x) - w(x)
    a, b = 0, 3
    e = 1e-3
    cmax = 50
    xo = racine.Dichotomie(d, a, b, e, cmax)[1]
    print(f"Xo = {xo:.4g} \t Yo ={w(xo):.4g}")
    
    xmin, xmax = -3, 3
    ymin, ymax = -10, 50
    graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    graphe.TraceAxes()
    graphe.TracePoints(X, Y)
    graphe.TraceFonc(f, xmin, xmax)
    graphe.TraceFonc(w, xmin, xmax, couleur='black')
    plt.show
# ==============================================================
def exo4():
    print("############# EXO 4 ##############")
    n, X, Y = disque.LireVec2('Data3.txt')
    
    # Chngmnt de Var
    X, Y = np.array(X), np.array(Y)
    V = gBETA*Y*Y
    sx = np.sin(X)
    U = gBETA*gBETA*(sx*(1 + sx) - 1)
    # plt.scatter(U, V)
    a, b = approx.Droite_MC(U, V)
    print(f"a = {a:.4g} \t b = {b:.4g}")
    def g(u):
        su = sin(u)
        return sqrt(a*gBETA*(su + su*su) + b/gBETA - a*gBETA)
    S = approx.EQM(X, Y, g)
    print(f"S = {S:.5g}")
    
    xmin, xmax = -3,3
    ymin, ymax = 0.4, 0.65
    graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    graphe.TraceAxes()
    graphe.TracePoints(X, Y)
    graphe.TraceFonc(g, xmin, xmax)
    plt.show
# ==============================================================
def exo5():
    pass
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================
