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
import syslin
import cbase
import graphe
import derive
import racine

# -----
gBETA = 0.1433
# ==============================================================
def main():
    exo3()

# ==============================================================
def exo1():
    print("############ EXO 1 ############")
    pr = 10 #nb de chiffre significatifs
    
    # Entrée des données
    a = 20*gBETA
    b = 40*gBETA
    A = [[0., 2., b, -1.],
         [2., 1., a, 1.],
         [-2, -1, a, 1.],
         [a, a, -1., 2.]]
    b = [-a, a, -2., b]
    A,b = np.array(A), np.array(b)
    # Résolution de sys lin
    x = syslin.Gauss_Jordan(A, b)
    print("\nLa solution est:")
    cbase.affiche_tableau([x], pr)
    
    # La matrice inverse
    invA = syslin.Inverse2(np.array(A))
    print("\nLa matrice inverse:")
    cbase.affiche_tableau(invA, pr)
    
    # Le determinant
    det = syslin.Determinant(A)
    print(f"\nLe determinant est:\n {det:.{pr}g}")
    
    # Le residu
    r = syslin.Residu(A, b, x)
    norme = syslin.Norme2(r)
    print("\nLe vecteur résidu:")
    cbase.affiche_tableau([r], pr)
    print(f"Sa norme: {norme:.{pr}g}")
# ==============================================================
def exo2():
    pass
# ==============================================================
def exo3():
    def getA(m):
        c = 50*gBETA/sqrt(m*m*m + 1)
        return [[1,c,-1],[2,c,1],[1,2,1]]
    
    def g(m):
        a = 50*gBETA
        b = [a,0,1]
        A = getA(m)
        x = syslin.Gauss_Jordan(A, b)
        return syslin.Norme2(x) - a
    
    def dg(m):
        return derive.dfc2(g, m, 0.001)
    lamd = racine.Dichotomie(dg, 0.7, 0.8, 0., 20)
    print(lamd)
    
    graphe.TraceFonc(g, 0, 5, npts=1001)
    plt.show()
# ==============================================================
def exo4():
    print("############ EXO 4 ############")
    pr = 10 #nb de chiffre significatifs
    
    for n in range(1,100,5):
        A,b = getAb(n)
        x = syslin.Gauss_Jordan(A,b)
        print(f"Pour n = {n} \t x =")
        cbase.affiche_tableau([x], pr)
    
    # print("\nPour n = 200")
    # A,b = getAb(n)
    # x = syslin.Gauss_Jordan(A, b)
    # norm = syslin.Norme2(x)
    # print(f"Norme de x = {norm:.{pr}}")
    
# ==============================================================
def exo5():
    pass
# ==============================================================
def getAb(n):
    A = np.zeros((n,n))
    b = np.zeros(n)
    for i in range(n):
        for j in range(n):
            A[i,j] = (i+j)*(i+j)
            b[i] += A[i,j]
    return A,b
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================
