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
import ctrl
import disque
import graphe
import syslin
import approx
import quadra
import racine
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
cs = 5         # Numbre of most significant digits to be printed
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo1()
    exo2()
    exo3()
    exo4()
    exo5()
# ==============================================================
def exo1():
    print("==================== EXO 1 =====================")

    A = [[4.,2,-3,1,0,2],[2,8,-2,0,-1,-1],[-1,-1,7,0,1,1],[0,-2,-1,5,-1,0],[-gBETA,1,1,2,6,3],[2,-gBETA,-1,1,-2,8]]
    b = [-6.,1,5,-7,10*gBETA,4]
    
    A,b = np.array(A),np.array(b)
    A_c = A.copy()
    
    x = syslin.Elim_Gauss0(A_c, b) # A_c sera modifiée
    print("x* =", end=' ')
    for comp in x:
        print(f"{comp:.{cs}g}", end=' ')
    print("\n")
    
    r = syslin.Residu(A, b, x)
    norme = syslin.Norme2(r)
    print(f"||R|| = {norme:.{cs}g}\n")
    
    det = syslin.Determinant(A)
    print(f"detA = {det:.{cs}g}\n")
    
    # invA = syslin.Inverse2(A) # Did not work well
    invA = np.linalg.inv(A)
    print(f"invA[2,2] = {invA[2,2]:.{cs}g}\t invA[4,3] = {invA[4,3]:.{cs}g}\n")
      
    A_3 = Gauss_modifie(A)
    cbase.affiche_tableau(A_3, cs, '  A à la 3ème étape  ')

    print("================================================\n")
# ==============================================================
def exo2():
    print("==================== EXO 2 =====================")
    
    n = 10 # nb inconnus
    # Mx = b  #   M tridiag
    L = np.zeros(n)
    D = np.zeros(n)
    U = np.zeros(n)
    b = np.zeros(n)
    for i in range(n):
        L[i] = g(1 + 10*i*gBETA)
        D[i] = 0.5*g(0.1 + sqrt(i)*gBETA)
        U[i] = g(0.3 + i*gBETA)
        b[i] = 0.5*i
    x = syslin.Gauss_Thomas(L, D, U, b)
    print(f"x1 = {x[0]:.{cs}g}\t x2 = {x[1]:.{cs}g}\t x3 = {x[2]:.{cs}g}\n")
    
    r = syslin.Residu_tridiag(L, D, U, b, x)
    norm = syslin.Norme2(r)
    print(f"||R|| = {norm:.{cs}}\n")
    
    print("================================================\n")
# ==============================================================
def exo3():
    print("==================== EXO 3 =====================")
    
    n, X, Y = disque.LireVec2('fonction.txt')
    X,Y = np.array(X), np.array(Y)
    
    U = np.log(X)
    V = Y - X*np.sin(X*X)
    c = approx.Poly_MC(U, V, 3)
    
    c0 = c[0]
    c1 = 2*c[1]/3
    c2 = 9*c[2]*log(10)*log(10)/16
    c3 = 64*c[3]/125
    print(f"c0 = {c0:.{cs}g}\t c1 = {c1:.{cs}g}\t c2 = {c2:.{cs}g}\t c3 = {c3:.{cs}g}\n")
    
    def ajust(x):
        return c0 + c1*log(x*sqrt(x)) + c2*log(x**(4/3),10)**2 + c3*log(x**(5/4))**3 + x*sin(x*x)
    
    eqm = approx.EQM(X, Y, ajust)
    print(f"EQM = {eqm:.{cs}}\n")
    
    xmin, xmax = 0.2,3.5
    ymin, ymax = -2.,8.
    graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    graphe.TraceAxes()
    graphe.TracePoints(X, Y)
    graphe.TraceFonc(ajust, xmin, xmax)
    
    # xmin, xmax = -2., 1.
    # ymin, ymax = -8, 8.
    # def Poly(x):
    #     return cbase.Horner(c, x)
    # graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    # graphe.TraceAxes()
    # graphe.TraceFonc(Poly, xmin, xmax, couleur='red')
    # plt.scatter(U,Y)
    print("================================================\n")
# ==============================================================
def exo4():
    print("==================== EXO 4 =====================")
    
    A = np.array([[4.,2,-3,1,0,2],[2,8,-2,0,-1,-1],[-1,-1,7,0,1,1],[0,-2,-1,5,-1,0],[-gBETA,1,1,2,6,3],[2,-gBETA,-1,1,-2,8]])
    xo = np.array([1.,0,0,0,0,0])
    iter_max = 1000
    eps = 1e-6
    q, n_iter = Stodola(A, iter_max, xo, eps)
    print(f"q = {q:.{cs}g}\t n_iter = {n_iter}\n")  
    # eig = np.linalg.eig(A)[0]
    # print(max(eig))
    
    print("================================================\n")
# ==============================================================
def exo5():
    print("==================== EXO 5 =====================")
    
    s = quadra.comp_Simpson2(sol_sys, -2., 3, 50)
    print(f"s = {s:.{cs}g}\n")    
    
    x = racine.Regula_Falsi(func, -2, 10, 1e-6, 150)[1]
    print(f"x = {x:.{cs}g}\n")
    
    # xmin, xmax = -2.,3.
    # ymin, ymax = -5.,5.
    # graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    # graphe.TraceAxes()
    # graphe.TraceFonc(sol_sys, xmin, xmax)
    
    print("================================================\n")
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def f(x):
    pass
# ==============================================================
def g(x): # pour exo2
    return sqrt(x) + exp(-x)
# ==============================================================
def Gauss_modifie(A): # pour exo1 Q5
    # Initiliazing Vars
    A = A.copy()
    n = len(A)
    
    # Triangularization of A
    for e in range(0,3): # arreter a la troisieme etape i = 2
        # This method assigns the pivot directly from the diagonal, with no row swapping
        pivot = A[e,e]
        if (pivot == 0):
            ctrl.erreur('syslin.Elim_Gauss0', 'pivot nul')
        # Eliminating terms under pivot
        for i in range(e+1,n):
            alpha = A[i,e]/pivot
            A[i,e] = 0.
            for j in range(e+1,n):
                A[i,j] = A[i,j] - alpha*A[e,j]
 
    return A[2], A[3] # La 3ème et 4ème ligne
# ==============================================================
def Stodola(A, iter_max, xo, eps):
    x = xo
    n_iter = 0
    while(n_iter < iter_max):
        y = A@x
        r = np.sqrt(sum(y*y))
        y /= r
        norm_d = syslin.Norme2(x - y)
        if (norm_d < eps):
            # return y, n_iter
            x = y
            break
        x = y
        n_iter += 1
    
    q = (x.T)@(A@x)/((x.T)@x)
    
    return q, n_iter
# ==============================================================
def sol_sys(t):
    A = np.array([[4.,2,-3,1,0,2],[2,8,-2,0,-1,-1],[-1,-1,7,0,1,1],[0,-2,-1,5,-1,0],[-gBETA,1,1,2,6,3],[2,-gBETA,-1,1,-2,8]])
    b = np.zeros(6)
    for i in range(6):
        b[i] = t - i/3
    x = syslin.Gauss_Jordan(A, b)
    
    return syslin.Norme2(x)
# ==============================================================
def func(x):
    return quadra.comp_Simpson2(sol_sys, -2, x, 50) - exp(-2*x) + 10*gBETA
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================