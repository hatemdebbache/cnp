"""
        Module : EDO
        Résolution des équation diff
"""
# ==========================================================================
import numpy as np
import syslin

"""
        Problème avec Conditions Aux Limites
        y'' = P*y' + Q*y + R
        
"""
# ==========================================================================
def DFC2_GT(foncPQR, a, b, Va, Vb, Nsub, condA, condB): # Ord 2
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
def DFC4_GT(foncPQR, a, b, Va, Vb, Nsub, condA, condB): # Ord 4
    x, y = DFC2_GT(foncPQR, a, b, Va, Vb, Nsub, condA, condB)
    x2, y2 = DFC2_GT(foncPQR, a, b, Va, Vb, Nsub*2, condA, condB)
    for i in range(Nsub+1):
        y[i] = y2[2*i] + (y2[2*i] - y[i])/3

    return x, y
# ==========================================================================

"""
        Problème de CAUCHY 1er ORD
        y' = f(t, y)

"""
# ==========================================================================
def Euler1(f, a, b, ya, Nsub, m): # Ord 1
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph.
    H = (b - a)/Nsub ; h = H/m
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
def Heun1(f, a, b, ya, Nsub, m): # Ord 2
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph
    H = (b - a)/Nsub ; h = H/m ; h2 = 0.5*h
    t = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    t[0] = a ; y[0] = ya
    v = ya
    for k in range(0,Nsub):
        u = t[k]
        for i in range(0,m):
            Po = f(u, v)
            P = f(u + h, v + h*Po)
            v = v + h2*(Po + P)
            u = u + h
        t[k+1] = a + (k+1)*H
        y[k+1] = v
    
    return 2*m*Nsub, t, y
# ==========================================================================
def PointMilieu1(f, a, b, ya, Nsub, m): # Ord 2
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph
    H = (b - a)/Nsub ; h = H/m ; h2 = 0.5*h
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
def RungeKutta1(f, a, b, ya, Nsub, m): # Ord 4
    H = (b - a)/Nsub ; h = H/m ; h2 = 0.5*h; h6 = h/6
    t = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    t[0] = a ; y[0] = ya
    v = ya
    for k in range(0,Nsub):
        u = t[k]
        for i in range(0,m):
            Po = f(u, v)
            Pm1 = f(u + h2, v + h2*Po)
            Pm2 = f(u + h2, v + h2*Pm1)
            P = f(u + h, v + h*Pm2)
            v = v + h6*(Po + 2.*(Pm1 + Pm2) + P)
            u = u + h
        t[k+1] = a + (k+1)*H
        y[k+1] = v
    
    return 4*m*Nsub, t, y
# ==========================================================================

"""
        Problème de CAUCHY 2ème ORD
        y', z' = f(t, y, z)
    
"""

# ==========================================================================
def Euler2(f, a, b, ya, za, Nsub, m): # Ord 1
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph.
    H = (b - a)/Nsub ; h = H/m
    t = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    z = np.zeros(Nsub+1)
    t[0] = a ; y[0] = ya ; z[0] = za
    v = ya ; w = za
    for k in range(0,Nsub):
        u = t[k]
        for i in range(0,m):
            Po, Qo = f(u, v, w)
            v = v + h*Po
            w = w + h*Qo
            u = u + h
        t[k+1] = a + (k+1)*H
        y[k+1] = v
        z[k+1] = w
    return m*Nsub, t, y, z
# ==========================================================================
def Heun2(f, a, b, ya, za, Nsub, m): # Ord 2
    # Nsub : int : graphique
    # m : int : nb de subdivision entre deux pas graph
    H = (b - a)/Nsub ; h = H/m ; h2 = 0.5*h
    t = np.zeros(Nsub+1)
    y = np.zeros(Nsub+1)
    z = np.zeros(Nsub+1)
    t[0] = a ; y[0] = ya ; z[0] = za
    v = ya ; w = za
    for k in range(0,Nsub):
        u = t[k]
        for i in range(0,m):
            Po, Qo = f(u,v,w)
            P, Q = f(u + h, v + h*Po, w + h*Qo)
            v = v + h2*(Po + P)
            w = w + h2*(Qo + Q)
            u = u + h
        t[k+1] = a + (k+1)*H
        y[k+1] = v
        z[k+1] = w
    return 2*m*Nsub, t, y, z
# ==========================================================================






# METHODES For 1st Ord ODE. IMPLEMENTED BY HATEM
# ==========================================================================
#
#                           O(h²)
#
# ==========================================================================
def Adams_Bashfort2(f, to, t1, tn, yo, y1):
    h = t1 - to ; h2 = 0.5*h
    t, y = [to, t1], [yo, y1]
    i = 1
    while(t[-1] < tn):
        k1 = f(t[i],y[i])
        k2 = f(t[i-1],y[i-1])
        t.append(to + (i+1)*h)
        y.append(y[i] + h2*(3.*k1 - k2))
        i += 1
    
    return 2*(i-1), np.array(t), np.array(y)
# ==========================================================================
def Euler_Mod(f, a, b, ya, Nsub, m):
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
            P1 = f(u + h2, v + h2*Po)
            v = v + h*P1
            u = u + h
        t[k+1] = a + (k+1)*H
        y[k+1] = v
    
    return 2*m*Nsub, t, y
# ==========================================================================
#
#                           O(h^3)
#
# ==========================================================================
def Adams_Bashfort3(f, to, t1, t2, tn, yo, y1, y2):
    h = t1 - to ; h12 = h/12
    t, y = [to, t1, t1], [yo, y1, y2]
    i = 2
    while(t[-1] < tn):
        k1 = f(t[i],y[i])
        k2 = f(t[i-1],y[i-1])
        k3 = f(t[i-2],y[i-2])
        t.append(to + (i+1)*h)
        y.append(y[i] + h12*(23.*k1 - 16.*k2 + 5.*k3))
        i += 1
    
    return 3*(i-2), np.array(t), np.array(y)
# ==========================================================================