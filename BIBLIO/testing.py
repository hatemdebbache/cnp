"""
    Liste des fonction pour tester les méthodes numérique

"""
import numpy as np
from math import exp, log, sin, cos, sqrt
import edo
gBETA = 0.1433


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    def simple(t, y):
        return -4*y + np.exp(-4*t)
    def sol(t):
        return t*np.exp(-4*t)
    # Setting environment
    plt.clf()
    a, b = 0., 5.
    num = 200
    m = 2
    # Analytic solution
    t = np.linspace(a, b, m*num*10)
    y = sol(t)
    plt.plot(t, y, color='blue', label = 'Analytic')
    
    # Euler1 solution
    c, t, y = edo.Euler1(simple, a, b, sol(a), num, m)
    plt.scatter(t, y, color='red', label='Euler-1', lw = 2.)
    
    # Heun1 solution
    c, t, y = edo.Heun1(simple, a, b, sol(a), num, m)
    plt.scatter(t, y, color='orange', label='Heun-1', lw = 2.)
    
    # PM solution
    c, t, y = edo.PointMilieu1(simple, a, b, sol(a), num, m)
    plt.scatter(t, y, color='yellow', label='PM', lw = 2.)
    
    # Euler_Mod solution
    c, t, y = edo.Euler_Mod(simple, a, b, sol(a), num, m)
    plt.scatter(t, y, color='green', label='Euler_Mod', lw = 2.)
    
    # RungeKutta solution
    c, t, y = edo.RungeKutta1(simple, a, b, sol(a), num, m)
    plt.scatter(t, y, color='black', label='RungeKutta', lw = 2.)
    
    # METHODS WITH MORE THAN INITIAL VALUE "multistep"
    h = (b - a)/(m*num)
    
    # Adams_Bashfort2 solution
    c, t, y = edo.Adams_Bashfort2(simple, a, a+h, b, sol(a), sol(a+h))
    plt.scatter(t, y, color='pink', label='AdamsBashfort-2', lw = 2.)
    
    # Adams_Bashfort2 solution
    c, t, y = edo.Adams_Bashfort3(simple, a, a+h, a+2.*h, b, sol(a), sol(a+h), sol(a+2.*h))
    plt.scatter(t, y, color='magenta', label='AdamsBashfort-3', lw = 2.)
    
    plt.legend(fontsize=16.)
    plt.show()
# ============================================================================
def f(x): # (beta*x)^4 - 11*(beta*x) + 8*beta
    bx = gBETA*x
    return (bx)*(bx)*(bx)*(bx) - 11*bx + 8*gBETA

def df(x):  # df(x)/dx
    bx = gBETA*x
    return gBETA*(4*(bx)*(bx)*(bx) - 11)

def d2f(x): # d2f(x)/dx2 === f''(x)
    bx = gBETA*x
    return 12*gBETA*gBETA*bx*bx

def d3f(x): # d3f(x)/dx3    
    return 24*gBETA*gBETA*gBETA*gBETA*x

def g(x): # f(x)/f'(x)
    bx = gBETA*x
    return (bx*bx*bx*x - 11*x + 8)/(4*bx*bx*bx - 11)

def w(x): # f(x) = 0 === w(x) = x
    bx = gBETA*x
    return (bx*bx*bx*x + 8)/11
# ============================================================================
def f1(x): # exp(beta*x) - 3*(beta*x)^2
    bx = gBETA*x
    return exp(bx) - 3*bx*bx

def df1(x):
    bx = gBETA*x
    return gBETA*(exp(bx) - 6*bx)

def d2f1(x):
    bx = gBETA*x
    return gBETA*gBETA*(exp(bx) - 6)

def d3f1(x):
    bx = gBETA*x
    return gBETA*gBETA*gBETA*exp(bx)

def g1(x):
    bx = gBETA*x
    return ( exp(bx) - 3*bx*bx )/( gBETA*(exp(bx) - 6*bx) )

def w1(x):
    bx = gBETA*x
    return sqrt(exp(bx)/(3*gBETA*gBETA))
# ============================================================================
def f2(x): # 2*sin(beta*x) - (beta*x)*log(beta*x)
    bx = gBETA*x
    return 2*sin(bx) - bx*log(bx)

def df2(x):
    bx = gBETA*x
    return gBETA*(2*cos(bx) - log(bx) - 1)

def d2f2(x):
    bx = gBETA*x
    return gBETA*(-2*gBETA*sin(bx) - 1/x)

def d3f2(x):
    bx = gBETA*x
    return gBETA*(-2*gBETA*gBETA*cos(bx) + 1/(x*x))

def g2(x):
    bx = gBETA*x
    return  ( 2*sin(bx) - bx*log(bx) )/( gBETA*(2*cos(bx) - log(bx) - 1) )

def w2(x):
    bx = gBETA*x
    # return 2*sin(bx)/(gBETA*log(bx))
    return ( 2*sin(bx) + bx*(1 - log(bx)) )/gBETA
# ============================================================================
