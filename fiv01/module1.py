# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020

@author: user
"""
# %matplotlib auto
# %matplotlib inline
from math import exp, pi, cos, sin, log
import numpy as np
import matplotlib.pyplot as plt
import sys
# -----
sys.path.insert(0, "./../BIBLIO/")
import graphe
import disque
import syslin

# -----
gBETA = 0.1433
# ==============================================================
def main():
    exo1()

# ==============================================================
def exo1():
    n,L,D,U,b = disque.LireTridiag('donne.txt')
    print(b)
    syslin.Gauss_Thomas2(L, D, U, b)
    print(b)
# ==============================================================
def exo2():
    b = [1, 5, -1, 10, 3., 5, 6, -4, 3, 4, 15]
    x = syslin.Gauss_Thomas_LDU(2, 4, 1, b)
    print(x)
# ==============================================================
def exo3():
    n,m,A = disque.LireMat('M1.txt')
    n,m,B = disque.LireMat('M2.txt')
    c = [1.2, 5.1, -3., 4.6]
    y = syslin.Elim_Gauss0(A, c)
    x = syslin.Elim_Gauss0(B, y)
    print(x)
# ==============================================================
def exo4():
    pass
# ==============================================================
def exo5():
    pass
# ==============================================================
def f(x):
    return (1+x*x)*exp(-gBETA*x)
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================
