# =============================================================================
def Dichotomie(f, a, b, cmax, eps):
    # CTRL DES ENTREES
    if (a > b):
        a, b = b, a
    v1 = f(a)
    v2 = f(b)
    c = 2
    if (v1*v2 > 0):
        print("Cette méthode ne peut être pas appliqué dans cet intervale")
        return
    # DEBUT
    d = abs(b - a)/2
    
    while (c < cmax and d > eps):
        v = f(a + d)
        if (v*v1 > 0):
            a = a + d
            v1 = v
        else:
            b = b - d
            v2 = v
        d /= 2
        c += 1
    
    if (d < eps):
        print("Opération terminé avec succes, le coût utilisé est: {}".format(c))
    else:
        print("Le coût max est consommé, l'éstimation de l'erreur est: {}".format(d))

    return b - d        
# =============================================================================        
def Regula_Falsi(f, a, b, cmax, eps):
    # CTRL DES ENTREES
    if (a > b):
        a, b = b, a
    v1 = f(a)
    v2 = f(b)
    c = 2
    if (v1*v2 > 0):
        print("Cette méthode ne peut être pas appliqué dans cet intervale")
        return
    # DEBUT
    e = -v1*(b - a)/(v2 - v1)
    
    while (c < cmax and abs(e) > eps):
        a = a + e
        v1 = f(a)
        c += 1
        e = -v1*(b - a)/(v2 - v1)
    
    if (abs(e) < eps):
        print("Opération terminé avec succes, le coût utilisé est: {}".format(c))
    else:
        print("Le coût max est consommé, l'éstimation de l'erreur est: {}".format(e))

    return a + e
# =============================================================================
def Newton_Raphson(f, a, b, cmax, eps):
    v = 0.
    
    return v
# =============================================================================
def Secante(f, a, b, cmax, eps):
    # CTRL DES ENTREES
    y1 = f(a)
    y2 = f(b)
    c = 2
    if (y1*y2 > 0):
        print("Cette méthode ne peut être pas appliqué dans cet intervale")
        return
    # DEBUT
    v = (a*y2 - b*y1)/(y2 - y1)
    e = v - b
    
    while (c < cmax and abs(e) > eps):
        b = v
        y2 = f(b)
        c += 1
        v = (a*y2 - b*y1)/(y2 - y1)
        e = v - b
        
    if (abs(e) < eps):
        print("Opération terminé avec succes, le coût utilisé est: {}".format(c))
    else:
        print("Le coût max est consommé, l'éstimation de l'erreur est: {}".format(e))

    return v
# =============================================================================

# FOR TESTING
if __name__ == "__main__":
    from math import sqrt
    
    def g(x):
        return - x*x + 5
    print(Regula_Falsi(g, 1, 6, 25, 1e-6))
    print(Dichotomie(g, -1, -6, 25, 1e-6))
    print(Secante(g, 1, 6, 25, 1e-6))
    