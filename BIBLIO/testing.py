"""
    Liste des fonction pour tester les méthodes numérique

"""
from math import exp, log, sin, cos, sqrt

gBETA = 0.1433

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
def f3(x):
    bx = gBETA*x
    return 

def df3(x):
    bx = gBETA*x
    return

def d2f3(x):
    bx = gBETA*x
    return 

def d3f3(x):
    bx = gBETA*x
    return

def g3(x):
    return 

def w3(x):
    bx = gBETA*x
    return 
# ============================================================================