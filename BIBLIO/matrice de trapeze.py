import numpy as np
from math import exp ,tan

def trapezcomp(f, a, b, n):
    """
    Composite trapezoidal function integration

    INPUTS:
    f:  the function to integrate
    a:  lower bound of integration
    b:  upper bound
    n:  number of panels to create between ``a`` and ``b``
    """

    # Initialization
    h = (b - a) / n
    x = a

    # Composite rule
    In = f(a)
    for k in range(1, n):
        x  = x + h
        In += 2*f(x)

    return (In + f(b))*h*0.5

def romberg(f, a, b, p):
    """
    Romberg integration

    INPUTS:
    f:  the function to integrate
    a:  lower bound of integration
    b:  upper bound
    p:  number of rows in the Romberg table
    """

    I = np.zeros((p, p))
    for k in range(0, p):
        # Composite trapezoidal rule for 2^k panels
        I[k, 0] = trapezcomp(f, a, b, 2**k)

        # Romberg recursive formula
        for j in range(0, k):
            I[k, j+1] = (4**(j+1) * I[k, j] - I[k-1, j]) / (4**(j+1) - 1)

        print(I[k,0:k+1])   # display intermediate results

    return I

def Rombergg(f, a, b, cmax, e): # If other Richardson methods are to be made, update R() and Tnew() and T[] only
    """
        Méthode d'intégration Romberg
        
        Utilise la méthode des trapèzes à pas variable, et en extrapolant ses résultats pour h tend vers 0. La méthode converge plus rapidement que la méthode composite des trapèzes

    Parameters
    ----------
    f : la fonction à intégrer, elle doit être de classe <= kmax
    
    a, b : intervale d'intégration 
    
    cmax : condition d'arrêt sur le coût de la méthode
    
    e : l'erreur désiré

    Returns
    -------
    c : le coût consommé par la méthode (0 si cmax est atteint)
    
    k : le nombre de pas prise pour estimer 
    
    x : la valeur de l'intégrale
    
    ee : estimation de l'erreur sur x

    """
    def Tnew(Told, a, b, Nsub):
        # Utilisé pour trouver la prochaine valeur de l'integral avec la methode des trapèzes
        # avec un pas double de Nsub, pour minimser le coût de calcule
        S = 0.
        H = (b - a)/Nsub
        for i in range(Nsub):
            S += f(a + (i + 0.5)*H)
    
        H = H/2.        
        
        return 0.5*Told + H*S
    
    def R(A, B, fac):
        # Le formule d'extrapolation de Richardson pour cette méthode élémentaire
        return A + (A - B)/fac
    
    # Initiliasing table T to store values of the integral at different steps
    # And evaluating the first step Nsub = 1, h = b - a
    T = [[0., 0.]]
    Nsub = 1
    T[0][0] = 0.5*(b - a)*(f(a) + f(b))
    T[0][1] = Tnew(T[0][0], a, b, Nsub)
    
    # Initiliasing returned variables
    c = 2
    k = 1
    x = T[0][1]
    ee = abs(T[0][1] - T[0][0])
    
    while (ee > e):    
        k += 1     
        # if (k > kmax): # maximum depth reached
        #     k = 0
        #     break
        
        # Exanding T to store the next level of extrapolation
        T.append([0.0, 0.0])
        
        for i in range(1, k - 1):
            T[i][1] = R(T[i - 1][1], T[i - 1][0], 2**(2*i) - 1)
        
        # Evaluating error on current step of extrapolation
        # Avoid the evaluation of the last column if not necessary
        ee = abs(T[k - 2][1] - T[k - 2][0])
        if (ee < e):
            x = T[k - 2][1]
            break
        
        T[k - 1][1] = R(T[k - 2][1], T[k - 2][0], 2**(2*(k - 1)) - 1)
        
        # Moving indexes, for minimal memory use
        for i in range(0, k):
            T[i][0] = T[i][1]
    
        # Checking current cost
        c += Nsub
        if (c > cmax):
            c = 0
            break
        
        # Evaluating the integral of the next step h = h/2
        Nsub *= 2
        T[0][1] = Tnew(T[0][0], a, b, Nsub)
        
        # Updating the returned value to the latest calculation
        x = T[k - 1][1]
        
    return c, k, T, ee

if __name__ == '__main__':
    def func(x):
        return exp(x)

    p_rows = 5
    I = romberg(func, 0, np.pi/2, p_rows)
    solution = I[p_rows-1, p_rows-1]
    print(solution)    
    print(Rombergg(func, 0, np.pi/2, 100, 1e-7)[2])