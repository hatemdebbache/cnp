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
import edo

# ----- Global Variables/Constants
gBETA = 0.1631 # PERSONAL CODE
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
    ya= 1+ gBETA
    a=0
    c,x,y =edo.PointMilieu1(f, a, 0.6, ya, 3, 1)
    cc,xx,yy=edo.RungeKutta1(f, a, 0.6, ya, 2, 1)
    ccc,xxx,yyy=edo.Euler1(f, a,0.6, ya, 6, 1)
    print(yyy[6])
    
    print(y[3])
    print(yy[2])
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
def f(x,y):
    return -2*y+x+4+2*gBETA
# ==============================================================
def g(x):
    pass
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================