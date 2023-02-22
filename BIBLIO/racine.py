"""
            Module: Racine
            Resolution des problèms non-linéaire de la forme f(x) = 0
            
            
            La signature générale des procédures
            (f, a, b, e, cmax) ===> PROCEDURE ===> (c, x, echap)
"""

import ctrl
from math import sqrt

inf = float('inf')
# =================================================================
#                       Methodes Fermées
# =================================================================
def Dichotomie(f, a, b, e, cmax):
    """

    Parameters
    ----------
    f : fonction de l'eq f(x) = 0
    
    [a, b]: le domaine ou il se trouve la solution x
    
    e : l'erreur désirée
    
    cmax : le coût maximal autorisé > 2

    Returns
    -------
    c : le coût consommé par la procedure
    x : l'approximation de la racine 
    e : estimation de l'erreur sur l'approximation

    """
    # CTRL DES ENTREES
    yn = f(a)
    yp = f(b)
    if (a == b):
        ctrl.erreur('racine.Dichotomie', '[a, b] machi m9bol')
    if (yn*yp > 0):
        ctrl.erreur('racine.Dichotomie', 'Cette méthode ne peut pas être appliquée dans cet intervale')
    if (e < 0):
        ctrl.erreur('racine.Dichotomie', 'Valeur d\'erreur non valide')        

    # DEBUT
    xn, xp = a, b
    if (yn > 0):
        xn, xp = xp, xn

    c = 2
    d = xp - xn
    while (True):         
        d = 0.5*d
        x = xp - d
        ee = abs(d)
        if (ee < e):    # Convergence
            break
        
        c = c + 1
        if (c > cmax): # Non convergence
            c = 0
            break

        y = f(x)
        if (y > 0):
            xp = x
            
    return c, x, ee
# =================================================================
def Regula_Falsi(f, a, b, e, cmax):
    """

    Parameters
    ----------
    f : fonction de l'eq f(x) = 0
    
    [a, b]: le domaine ou il se trouve la solution x
    
    e : l'erreur désirée
    
    cmax : le coût maximal autorisé

    Returns
    -------
    c : le coût consommé par la procedure
    
    x : l'approximation de la racine 
    
    e : estimation de l'erreur sur l'approximation

    """
    # CTRL DES ENTREES
    yn = f(a)
    yp = f(b)
    if (a == b):
        ctrl.erreur('racine.Regula_Falsi', '[a, b] machi m9bol')
    if (yn*yp > 0):
        ctrl.erreur('racine.Regula_Falsi', 'Cette méthode ne peut pas être appliquée dans cet intervale')
    if (e < 0):
        ctrl.erreur('racine.Regula_Falsi', 'Valeur d\'erreur non valide')        

    # DEBUT
    xn, xp = a, b
    if (yn > 0):
        xn, xp = xp, xn
        yn, yp = yp, yn
    c = 2
    ee = inf
    x = xn
    
    while (True):        
        anc = x
        x = xn - yn*(xp - xn)/(yp - yn)
        ee = abs(x - anc)
        if (ee < e):    # Convergence
            break

        c = c + 1
        if (c > cmax): # Non convergence
            c = 0
            break
        
        y = f(x)
        if (y < 0):
            xn, yn = x, y
        else:
            xp, yp = x, y
            
    return c, x, ee
# =================================================================
#                       Methodes ouvertes 
# =================================================================
def Newton_Raphson(g, a, b, e, cmax, xo):
    """

    Parameters
    ----------
    g : La fonction f(x)/f'(x)
    [a, b] : le domaine ou la solution existe (utilisé pour détecter la divergence)
    e : l'erreur désiré
    cmax : le coût maximal autorisé
    xo : point de départ de la suite

    Returns
    -------
    c : le cout utilisé
        0 dans le cas de la non-convergence (cmax est consommé)
        <0 dans le cas de la divergence
    x : la solution approchée
        inf dans le cas de divergence
    ee : estimation de l'erreur sur la solution
        inf dans la cas de divergence

    """
    # CTRL
    if (a == b):
        ctrl.erreur('racine.Newton_Raphson', '[a, b] machi m9bol')
    if ((xo - a)*(xo - b) > 0):
        ctrl.erreur('racine.Newton_Raphson', 'xo machi m9bol')    
    if (e < 0):
        ctrl.erreur('racine.Newton_Raphson', 'Valeur d\'erreur non valide')        

    # DEBUT 
    c = 0
    x = xo
    ee = inf
    
    while (ee > e):
        c = c + 1
        if (c > cmax): # Non convergence
            c = 0
            break
        d = g(x)
        x = x - d
        if ((x - a)*(x - b) > 0):     # Divergence
            c, x, ee = -c, inf, inf
            break
        ee = abs(d)
    
    return c, x, ee
# =================================================================
def Secante(f, a, b, e, cmax, xo, x1):
    """

    Parameters
    ----------
    f : la fonction qu'on désire trouver une racine x
    [a, b]: le domaine qui délimte la racine
    e : l'erreur désirée
    cmax : le coût maximale autorisé
    x0, x1: pts de départ nécessaire pour cette méthode

    Returns
    -------
    c : le coût consommé
    x : la valeur approchée de la racine
    e : estimation de l'erreur sur l'approximation

    """
    # CTRL DES ENTREES
    if (a == b):
        ctrl.erreur('racine.Secante', '[a, b] machi m9bol')
    if (e < 0):
        ctrl.erreur('racine.Secante', 'Valeur d\'erreur non valide')        
    
    # DEBUT
    ee = inf
    yo = f(xo)
    y1 = f(x1)
    c = 2
    x = x1 # Initialisation de x
    d = x1 - xo
    while (True):
        d = -y1*d/(y1 - yo)
        x = x + d        
        if ((x - a)*(x - b) > 0): # Divergence
            c, x, ee = -c, inf, inf
            break
        
        ee = abs(d) 
        if (ee < e):  # Convergence
            break
        
        c = c + 1
        if (c > cmax):  # Non convergence
            c = 0
            break

        yo = y1
        y1 = f(x)
        
    return c, x, ee
# =================================================================
def Point_Fixe(w, a, b, e, cmax, xo):
    # CTRL
    if (a == b):
        ctrl.erreur('racine.Point_Fixe', '[a, b] machi m9bol')
    if ((xo - a)*(xo - b) > 0):
        ctrl.erreur('racine.Point_Fixe', 'xo machi m9bol')    
    if (e < 0):
        ctrl.erreur('racine.Point_Fixe', 'Valeur d\'erreur non valide')        

    # DEBUT 
    c = 0
    x = xo
    ee = inf
    
    while (ee > e):
        c = c + 1
        if (c > cmax): # Non convergence
            c = 0
            break
        
        anc = x
        x = w(x)
        
        if ((x - a)*(x - b) > 0):     # Divergence
            c, x, ee = -c, inf, inf
            break
        
        ee = abs(x - anc)
    
    return c, x, ee
# =================================================================
def racine_avec_poids(f, a, b, e, cmax):
    """ 
        New method that solves g(x) = f(x)*exp(m*(x - a)/h)
        The roots of g(x) are the same as f(x)
    """
    # CTRL DES ENTREES
    yn = f(a)
    yp = f(b)
    if (a == b):
        ctrl.erreur('racine.New', '[a, b] machi m9bol')
    if (yn*yp > 0):
        ctrl.erreur('racine.New', 'Cette méthode ne peut pas être appliquée dans cet intervale')
    if (e < 0):
        ctrl.erreur('racine.New', 'Valeur d\'erreur non valide')        

    # DEBUT
    x1, x2 = a, b
    y1, y2 = f(a), f(b)
    h = 0.5*(x2 - x1)
    x = 0.5*(x2 + x1)
    y = f(x)
    c = 3
    
    while (True):
        xnew = x - h*y/sqrt(y*y - y1*y2)
        ee = abs(xnew - x)
        if (ee < e):
            break
        c += 2
        if (c > cmax):
            c = 0
            break
        x1, x2 = x, xnew
        y1, y2 = y, f(xnew)
        h = 0.5*(xnew - x)
        x = 0.5*(xnew + x)
        y = f(x)
        
    return c, x, ee