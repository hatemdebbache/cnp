"""
        Module Approx
        Méthodes pour approximer nuage de pts avec des modèles mathématiques
"""
import numpy as np
from ctrl import erreur

# =============================================================================
def Droite_MC(X, Y, kmin = 0, kmax = None):
    """
    Least Squared Regression of data points (X, Y) in form of linear equation y = a*x + b.
    If only a portion of the data is used, indexes (kmin, kmax) should be provided.
    
    Parameters
    ----------
    X : array_like
        First axis coordinates of data.
    Y : array_like
        Second axis coordinates of data.
    kmin : int, optional
        Index of first element. The default is 0.
    kmax : int, optional
        (Included) Index of last element. The default is len(X).
        
    Returns
    -------
    a, b : float
        Coefficients of linear regression function (Y = a*X+b).

    """
    # Input Control
    n = len(X) 
    if n != len(Y):
        erreur('Droite_MC', 'Arrays X, Y must be of same dimension')
    if kmax == None:
        kmax = n-1
    elif kmax >= n:
        erreur('Droite_MC', 'Inavlid interval: kmax out of range')
    if kmin > kmax:
        erreur('Droite_MC', 'Invalid interval: kmax must be > kmin')
    
    # Calculation
    Xm, Ym = float(sum(X[kmin:kmax+1]))/(kmax+1-kmin), float(sum(Y[kmin:kmax+1]))/(kmax+1-kmin) # Cast to float is forced to avoid euclidian division
    Cov, VarX = 0.0,0.0
    for x, y in zip(X[kmin:kmax+1], Y[kmin:kmax+1]):
        Dx = x - Xm
        Cov += Dx*(y - Ym)
        VarX += Dx*Dx
        
    a = Cov/VarX
    b = Ym - a*Xm
    
    return a, b
# =============================================================================
def Poly_MC(X, Y, deg, kmin=0, kmax=None):
    """
    Least Squared Regression of data points (X,Y) in form of polynomial of degree 'deg'
        P(x) = c0 + c1*x + c2*x^2 + ... + cm*x^m  / m = deg.
    If only a portion of the data is used, indexes (kmin, kmax) should be provided.
    
    Parameters
    ----------
    X : array_like
        First axis coordinates of data.
    Y : array_like
        Second axis coordinates of data.
    deg : int
        Degree of regression polynomial.
    kmin : int, optional
        Index of first element. The default is 0.
    kmax : int, optional
        (Included) Index of last element. The default is len(X).

    Returns
    -------
    C : numpy.ndarray
        Coeffients of resulting polynomial.

    """
    # Input Control
    n = len(X) 
    if n != len(Y):
        erreur('Poly_MC', 'Arrays X, Y must be of same dimension')
    if type(deg) != int:
        erreur('Poly_MC', 'Polynomial\'s degree must be integer')
    if deg <= 0:
        erreur('Poly_MC', 'Polynomial\'s degree must be > 0')
    if kmax == None:
        kmax = n-1
    elif kmax >= n:
        erreur('Poly_MC', 'Inavlid interval: kmax out of range')
    if kmin > kmax:
        erreur('Poly_MC', 'Invalid interval: kmax must be > kmin')
    
    # Initialising Vars
    npts = deg+1
    A = np.zeros((npts, npts))
    b = np.zeros(npts)
    x_copy = np.array(X[kmin:kmax+1].copy()) # Converted to ndarray to profit from numpy.mul
    y_copy = np.array(Y[kmin:kmax+1].copy()) # Copy instead of modifying X, Y
    
    # Calculating
    b[0] = sum(y_copy)
    A[0,0] = n
    # Filling the left-lower L of A
    for i in range(1,npts):
        b[i] = sum(x_copy*y_copy)
        A[i,0] = sum(x_copy)
        x_copy *= X[kmin:kmax+1]
    for j in range(1,npts):
        A[npts-1,j] = sum(x_copy)
        x_copy *= X[kmin:kmax+1]
        
    # Filling the rest of the system (anit-diaognally)
    for j in range(1,npts):
        for i in range(npts-1):
            A[i,j] = A[i+1,j-1]
    
    # Evaluating the array of coeffs
    c = np.linalg.solve(A,b)
    
    return c
# =============================================================================
def EQM(X, Y, f=None, Y_approch=None, kmin = 0, kmax = None): 
    """
    Evaluate the Mean Squared Error (écart quadratique moyen) between data points (X, Y) and function f(x) (or data points (X, Y_approch)).

    Parameters
    ----------
    X : array_like
        first axis coordinates of data.
    Y : array_like
        second axis coordinates of data.
    f : function, optional
        The function of mathematical model used to approach data points. The default is None.
    Y_approch : array_like, optional
        The second axis coordinates of approching model points. The default is None.    
    kmin : int, optional
        Index of first element. The default is 0.
    kmax : int, optional
        (Included) Index of last element. The default is len(X).

    Note: at least one of the two arguments (f, Y_approch) must be provided. If both arguments are provided, f will overwrite Y_approch.
    
    Returns
    -------
    eqm : float
        Numerical value of MSE.

    """
    # Input Control
    n = len(X) 
    if n != len(Y):
        erreur('EQM', 'Arrays X, Y must be of same dimension')
    if kmax == None:
        kmax = n-1
    elif kmax >= n:
        erreur('EQM', 'Inavlid interval: kmax out of range')
    if kmin > kmax:
        erreur('EQM', 'Invalid interval: kmax must be > kmin')
    if f == None:
        if Y_approch.all == None:
            erreur('EQM', 'At least one of two arguments (f, Y_approch) must be provided')
        elif len(Y_approch) != (kmax+1-kmin):
            erreur('EQM', 'Y_approch must be of same dimension as the considered portion of (X, Y)')
    else:
        Y_approch = np.zeros(kmax+1-kmin)
        for i in range(kmin, kmax+1):
            Y_approch[i-kmin] = f(X[i])
    
    eqm = Y[kmin:kmax+1]-Y_approch
    eqm = np.sqrt(sum(eqm*eqm)/(kmax+1-kmin))
    
    return eqm
# =============================================================================