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
import disque
import syslin

# -----
gBETA = 0.1433
# ==============================================================
def main():
    exo1()

# ==============================================================
def exo1():
    n, m, A, b = disque.LireSys('system.txt')
    x = syslin.Elim_Gauss0(A, b)
    
    print(x)
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
def f(x):
    return (1+x*x)*exp(-gBETA*x)
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================
