"""
        Module : EDO
        Résolution des équation diff
"""
# ==========================================================================
import numpy as np
import ctrl
import syslin
# ==========================================================================
def DFC2_GT(foncPQR, a, b, Va, Vb, Nsub, condA, condB):
    x = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    m = Nsub - 1
    L = np.zeros(m)
    D = np.zeros(m)
    U = np.zeros(m)
    w = np.zeros(m)
    # Init
    h = (b - a)/Nsub; hh = h*h; h2 = 0.5*h
    # Remplir system
    for k in range(0,m):
        t = a + (k + 1)*h
        x[k+1] = t
        P, Q, R = foncPQR(t)
        D[k] = -2. - hh*Q
        t = h2*P
        L[k] = 1. + t
        U[k] = 1. - t
        w[k] = hh*R
    # Correction a gauche
    if (condA == True):
        w[0] = w[0] - L[0]*Va
    else:
        t = L[0]/3.
        w[0] = w[0] + 2.*h*Va*t
        D[0] = D[0] + 4.*t
        U[0] = U[0] - t
    # Correction a droite
    if (condB == True):
        w[m-1] = w[m-1] - U[m-1]*Vb
    else:
        t = U[m-1]/3.
        w[m-1] = w[m-1] - 2.*h*Vb*t
        L[m-1] = L[m-1] - t
        D[m-1] = D[m-1] + 4.*t
    # Resolution du sys
    syslin.Gauss_Thomas2(L, D, U, w)
    # Remplir y
    for k in range(1,Nsub):
        y[k] = w[k-1]
    # Limite a gauche
    x[0] = a
    if (condA == True):
        y[0] = Va    
    else:
        y[0] = (4.*y[1] - y[2] - 2.*h*Va)/3.
    # Limite a droite
    x[Nsub] = b
    if (condB == True):
        y[Nsub] = Vb
    else:
        y[Nsub] = (4.*y[Nsub-1] - y[Nsub-2] + 2.*h*Vb)/3.
    
    return x, y
# ==========================================================================
def DFC4_GT(foncPQR, a, b, Va, Vb, Nsub, condA, condB):
    x, y = DFC2_GT(foncPQR, a, b, Va, Vb, Nsub, condA, condB)
    x2, y2 = DFC2_GT(foncPQR, a, b, Va, Vb, Nsub*2, condA, condB)
    for i in range(Nsub+1):
        y[i] = y2[2*i] + (y2[2*i] - y[i])/3

    return x, y
# ==========================================================================

"""
                   CAUCHY 1er ORD
        y' = f(t, y)

"""
# ==========================================================================
def Euler1(f, a, b, ya, Nsub, m):
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph
    H = (b - a)/Nsub
    h = H/m
    
    t = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    
    t[0] = a ; y[0] = ya
    v = ya
    for k in range(0,Nsub):
        u = t[k]
        for i in range(0,m):
            v = v + h*f(u, v)
            u = u + h

        t[k+1] = a + (k+1)*H
        y[k+1] = v
    
    return m*Nsub, t, y
# ==========================================================================
def Heun1(f, a, b, ya, Nsub, m):
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph
    H = (b - a)/Nsub
    h = H/m ; h2 = 0.5*h
    
    t = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    
    t[0] = a ; y[0] = ya
    v = ya
    for k in range(0,Nsub):
        u = t[k]
        for i in range(0,m):
            Po = f(u, v)
            P1 = f(u + h, v + h*Po)
            v = v + h2*(Po + P1)
            u = u + h
        t[k+1] = a + (k+1)*H
        y[k+1] = v
    
    return 2*m*Nsub, t, y
# ==========================================================================
def PointMilieu1(f, a, b, ya, Nsub, m):
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph
    H = (b - a)/Nsub
    h = H/m ; h2 = 0.5*h
    
    t = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    
    t[0] = a ; y[0] = ya
    v = ya
    for k in range(0,Nsub):
        u = t[k]
        for i in range(0,m):
            Po = f(u, v)
            Pm = f(u + h2, v + h2*Po)
            v = v + h*Pm
            u = u + h
        t[k+1] = a + (k+1)*H
        y[k+1] = v
    
    return 2*m*Nsub, t, y
# ==========================================================================
# ==========================================================================
# ==========================================================================