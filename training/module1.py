# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020

@author: user
"""
# %matplotlib auto
# %matplotlib inline
from math import exp, pi, sin, log, cos
import numpy as np
import matplotlib.pyplot as plt
import sys
# -----
sys.path.insert(0, "./../BIBLIO/")
import graphe
import disque
from racine import *
from derive import * 
from quadra import *


# -----
gBETA = 0.1433
# ==============================================================
def main():
    exo2()
# ==============================================================
def exo1():
    print(Dichotomie(y, 0, 5, 1e-9, 20))
    print(Newton_Raphson(yy, 0, 5, 1e-9, 20, 2))
# ==============================================================
def exo2():
    v = 4000*gBETA*cos(800)*sin(800)
    h = 0.0001  # 1e-3  ==> 1e-15 2e-k = 1e-12
    while(h > 1e-15):
        w = dfc2(j, 800, h)
        h = 0.5*h
        print(abs(w - v))
# ==============================================================
def exo3():
    pass
# ==============================================================
def exo4():
    pass
# ==============================================================
def f(x):
    return log(x + 1) + 1 - x
def g(x):
    return -(x + 1)*(log(x + 1) + 1 - x)/x
# ==============================================================
def y(x):
    bx = gBETA*x
    return bx**5 - 13*bx**3 + 5*gBETA

def yy(x):
    bx = gBETA*x
    return (bx**5 - 13*bx**3 + 5*gBETA)/(5*gBETA*bx**4 - 39*gBETA*bx**2)
# ==============================================================
def j(x):
    return 2000*gBETA*(sin(x)**2)
# ==============================================================
def o(x):
    return exp(x*sin(x)) - (sin(x)/x)**3
# ==============================================================
if (__name__ == "__main__"):
    
    # main()
    gBETA = 0.1071
    X = [.25 + i for i in range(1,6)]
    y = [3, 2-2*gBETA, 2 - 4*gBETA/3 ,2, 1]
    print(y[3] - 2*y[2] + y[1])
    
# ==============================================================
