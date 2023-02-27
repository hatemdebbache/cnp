"""
    Module syslin
    Resolution des systems lineaires Ax=b
"""
import numpy as np
import ctrl


# =============================================================================
#
#                  SYSTEMES DENSES            
# 
# =============================================================================
def Elim_Gauss0(A, b):
    """
    Compute the solution of linear system Ax = b with Gauss-Jordan method *without* pivoting.

    Parameters
    ----------
    A : array_like
        Coefficient matrix.
    b : array_like
        Right hand side of the linear system.

    Returns
    -------
    x : ndarray
        Approximate solution.

    See Also:
        Gauss_Jordan
    
    """
    # Input Control
    n = len(b)
    if A.shape != (n,n):
        ctrl.erreur('syslin.Elim_Gauss0', f'A must be a square matrix of same dimension as b (expected ({n},{n}), got {A.shape}).')
    
    # Initiliazing Vars
    x = np.zeros(n)
    
    # Triangularization of A
    for e in range(0,n-1):
        # This method assigns the pivot directly from the diagonal, with no row swapping
        pivot = A[e,e]
        if (pivot == 0):
            ctrl.erreur('syslin.Elim_Gauss0', 'pivot nul')
        # Eliminating terms under pivot
        for i in range(e+1,n):
            alpha = A[i,e]/pivot
            A[i,e] = 0.
            for j in range(e+1,n):
                A[i,j] = A[i,j] - alpha*A[e,j]
            b[i] = b[i] - alpha*b[e]
    
    # Computing x with Ascendance
    x[n-1] = b[n-1]/A[n-1,n-1]
    for i in range(n-2,-1,-1):
        s = b[i]
        for j in range(i+1,n):
            s = s - A[i,j]*x[j]
        x[i] = s/A[i,i]
    
    return x
# ============================================================================
def Gauss_Jordan(A, b):
    """
    Compute the solution of linear system Ax = b with Gauss-Jordan's Elimination with partial pivoting.
    
    Parameters
    ----------
    A : array_like
        Coefficients Matrix.
    b : array_like
        Right Hand Side of the linear system.
    
    Returns
    -------
    x : ndarray
        Approximate solution of the linear system. 
    
    See Also:
        Gauss_Elim0
        
    """
    # Input Control
    A,b = check_ndarray(A,b)
    n = len(b)
    if A.shape != (n,n):
        ctrl.erreur('syslin.Gauss', 'A must be a square matrix of same order as b (expected ({n},{n}), got {A.shape}).')
    
    # Initialize Vars
    x = np.zeros(n)
    A = A.copy() # Copy to avoid altering
    b = b.copy()

    for e in range(n-1):
        Pivotation(A,n,e,b) # Pivotation before real calculation
        pivot = A[e,e]
        # if (pivot < 1e-16): ctrl.erreur('syslin.Gaus_Jordan', 'Singular Matrix')
        # Gauss-Jordan's elimination
        for i in range(e+1,n):
            coeff = A[i,e]/pivot
            b[i] -= coeff*b[e]
            A[i] -= coeff*A[e]
    
    # Compute the soltion by ascending
    x[n-1] = b[n-1]/A[n-1,n-1]
    for i in range(n-2,-1,-1):
        summ = b[i]
        for j in range(i+1,n):
            summ -= A[i,j]*x[j]
        x[i] = summ/A[i,i]

    return x
# ============================================================================
#
#                   SYSTEM TRIDIAG
#
# ============================================================================
def Gauss_Thomas(L, D, U, b):
    """
    Compute the solution of linear system Ax = b where A is a tridiagonal matrix.
    (Transform the tridiagonal matrix to upper bidiagonal matrix and compute each unknown recursively starting with the last).
    
    Parameters
    ----------
    L : array_like
        Lower Diagonal Coefficients.
    D : array_like
        Diagonal Coefficients.
    U : array_like
        Upper Diagonal Coefficients.
    b : array_like
        Right Hand Side of the linear system.

    Returns
    -------
    x : ndarray
        Approximate solution of the linear system.

    """
    n = len(b)
    
    # Initializing Vars
    x = np.zeros(n)
    c = np.zeros(n)
    
    # Computing Bidiagonal Matrix Coefficients (stored in x and c)
    p = D[0]
    x[0] = b[0]/p
    for i in range(1,n):
        m = i - 1
        c[m] = U[m]/p
        p = D[i] - L[i]*c[m]
        if (p == 0): ctrl.erreur('syslin.Gauss_Thomas', 'pivot nul')
        x[i] = (b[i] - L[i]*x[m])/p
        
    # Computing solution with ascendance
    for i in range(n-2,-1,-1):
        x[i] = x[i] - c[i]*x[i+1]
    
    return x
# ============================================================================
def Gauss_Thomas_D(L, D, U, b):
    """
    Compute the solution of linear system Ax = b where A is a tridiagonal matrix.

    Parameters
    ----------
    L : array_like
        Lower Diagonal Coefficients.
    D : float
        Digonal (all coefficients are equal to D).
    U : array_like
        Upper Diagonal Coefficients.
    b : array_like
        Right Hand Side of the linear system.

    Returns
    -------
    x : ndarray
        Approximate solution of the linear system.

    """
    n = len(b)
    
    # Initializing Vars
    x = np.zeros(n)
    c = np.zeros(n)
    
    # Computing Bidiagonal Matrix Coefficients (stored in x and c)
    p = D
    x[0] = b[0]/p
    for i in range(1,n):
        m = i - 1
        c[m] = U[m]/p
        p = D - L[i]*c[m]
        if (p == 0): ctrl.erreur('syslin.Gauss_Thomas_D', 'pivot nul')
        x[i] = (b[i] - L[i]*x[m])/p
        
    # Computing solution with ascendance
    for i in range(n-2,-1,-1):
        x[i] = x[i] - c[i]*x[i+1]
    
    return x
# ============================================================================
def Gauss_Thomas_LU(L, D, U, b):
    """
    Compute the solution of linear system Ax = b where A is a tridiagonal matrix.

    Parameters
    ----------
    L : float
        Lower Diagonal (All coefficients are equal to L).
    D : array_like
        Diagonal Coefficients.
    U : float
        Upper Diagonal (All coefficients are equal to U).
    b : array_like
        Right Hand Side of the linear system.

    Returns
    -------
    x : ndarray
        Approximate solution of the linear system.

    """
    n = len(b)
    
    # Initializing Vars
    x = np.zeros(n)
    c = np.zeros(n)
    
    # Computing Bidiagonal Matrix Coefficients (stored in x and c)
    p = D[0]
    x[0] = b[0]/p
    for i in range(1,n):
        m = i - 1
        c[m] = U/p
        p = D[i] - L*c[m]
        if (p == 0): ctrl.erreur('syslin.Gauss_Thomas_LU', 'pivot nul')
        x[i] = (b[i] - L*x[m])/p
        
    # Computing solution with ascendance
    for i in range(n-2,-1,-1):
        x[i] = x[i] - c[i]*x[i+1]
    
    return x
# ============================================================================
def Gauss_Thomas_LDU(L, D, U, b):
    """
    Compute the solution of linear system Ax = b where A is a tridiagonal matrix.

    Parameters
    ----------
    L : float
        Lower Diagonal (All coefficients are equal to L).
    D : float
        Diagonal (All coefficients are equal to D).
    U : float
        Upper Diagonal (All coefficients are equal to U).
    b : array_like
        Right Hand Side of the linear system.

    Returns
    -------
    x : ndarray
        Approximate solution of the linear system.

    """
    n = len(b)
    
    # Initializing Vars
    x = np.zeros(n)
    c = np.zeros(n)
    
    # Computing Bidiagonal Matrix Coefficients (stored in x and c)
    p = D
    x[0] = b[0]/p
    for i in range(1,n):
        m = i - 1
        c[m] = U/p
        p = D - L*c[m]
        if (p == 0): ctrl.erreur('syslin.Gauss_Thomas_LDU', 'pivot nul')
        x[i] = (b[i] - L*x[m])/p
        
    # Computing solution with ascendance
    for i in range(n-2,-1,-1):
        x[i] = x[i] - c[i]*x[i+1]
    
    return x
# ============================================================================
def Gauss_Thomas2(L, D, U, b):
    """
    Compute the solution of linear system Ax = b where A is a tridiagonal matrix.
    (Transform the tridiagonal matrix to upper bidiagonal matrix and compute each unknown recursively starting with the last).
    -**Enhaced**-
    
    Parameters
    ----------
    L : array_like
        Lower Diagonal Coefficients.
    D : array_like
        Diagonal Coefficients.
    U : array_like
        Upper Diagonal Coefficients.
    b : array_like
        Right Hand Side of the linear system.

    Returns
    -------
    x : ndarray
        Approximate solution of the linear system.
    
    See also:
        Gauss_Thomas
        
    """
    n = len(b)
    
    # Initializing Vars
    c = np.zeros(n)
    
    # Computing Bidiagonal Matrix Coefficients (stored in x and c)
    p = D[0]
    b[0] = b[0]/p
    for i in range(1,n):
        m = i - 1
        c[m] = U[m]/p
        p = D[i] - L[i]*c[m]
        if (p == 0): ctrl.erreur('syslin.Gauss_Thomas', 'pivot nul')
        b[i] = (b[i] - L[i]*b[m])/p
        
    # Computing solution with ascendance
    for i in range(n-2,-1,-1):
        b[i] = b[i] - c[i]*b[i+1]
# ============================================================================
#
#                   Supplementary Functions
#
# ============================================================================
def Determinant(A):
    """
    Compute the determinant of a matrix.

    Parameters
    ----------
    A : array_like
        *Square* matrix of which the determinant is computed.

    Returns
    -------
    d : float
        The Determinant of A.

    """
    # Input Control
    A = check_ndarray(A)
    try:
        n, m = A.shape
    except ValueError:
        ctrl.erreur('syslin.Determinant', 'A must be an array of 2-axis')
    if n != m:
        ctrl.erreur('syslin.Determinant', 'Non square matrices cannot be inverted')
    
    # Initializing Vars
    A = A.copy() # Copying A to avoid altering it
    d = 1.
    
    # Triangularisation of A
    for e in range(0,n-1):
        altered = Pivotation(A,n,e)
        pivot = A[e,e] 
        # if (pivot < 1e-15): d = 0.; break # Break loop for null pivot
        
        d *= pivot if not altered else -pivot # Diagonal product
        
        for i in range(e+1,n):
            alpha = A[i,e]/pivot
            A[i,e] = 0.
            for j in range(e+1,n): # Elements under digonal are ignored, won't affect determinant
                A[i,j] = A[i,j] - alpha*A[e,j]        
    d *= A[n-1,n-1] # Last digonal element    
    
    return d
# =============================================================================
def Inverse(A):
    """
    Compute the (multiplicative) inverse matrix of A

    Parameters
    ----------
    A : array_like
        Martix to be inversed.

    Returns
    -------
    invA : ndarray
        The inverse matrix of A.

    """
    # Input Control
    A = check_ndarray(A)
    try:
        n, m = A.shape
    except ValueError:
        ctrl.erreur('syslin.Inverse', 'A must be an array of 2-axis')
    if n != m:
        ctrl.erreur('syslin.Inverse', 'Non square matrices cannot be inverted')
    
    # Initializing var
    A = A.copy() # Copying A to avoid altering it
    invA = np.identity(n)
    
    # Transforming A to triangular (sup) matrix
    for e in range(0,n-1):
        pivot = A[e,e]        
        if (pivot == 0): ctrl.erreur('syslin.Inverse', 'pivot nul durant triangularisation')        
        for i in range(e+1,n):
            alpha = A[i,e]/pivot
            A[i,e] = 0.
            for j in range(e+1,n):
                A[i,j] = A[i,j] - alpha*A[e,j]
            for j in range(0,e+1):
                invA[i,j] = invA[i,j] - alpha*invA[e,j]
    
    # Transforming A to diagonal matrix
    for e in range(1, n):
        pivot = A[e,e]
        if (pivot == 0):
            ctrl.erreur('syslin.Inverse', 'pivot nul durant diagonalisation')
        for i in range(0, e):
            alpha = A[i,e]/pivot
            A[i,e] = 0.
            for j in range(e+1,n):
                A[i,j] = A[i,j] - alpha*A[e,j]
            for j in range(0,n):
                invA[i,j] = invA[i,j] - alpha*invA[e,j]
    
    # Transforming A to identity matrix
    for i in range(n):
        for j in range(n):
            invA[i,j] = invA[i,j]/A[i,i]
    
    return invA # Matrice inverse de A
# =============================================================================
def Residu(A, b, x):
    """
    Compute the residual of an estimation x for the solution of Ax = b 

    Parameters
    ----------
    A : array_like
        Coefficient matrix.
    b : array_like
        Right hand side of the linear system.
    x : array_like
        The estimation of the exact solution.

    Returns
    -------
    r : ndarray
        The *noramlized* estimation residual (Normalized = end result divided by vector length).

    """
    # Input Control
    A,b,x = check_ndarray(A,b,x)
    try:
        n, m = A.shape
    except ValueError:
        ctrl.erreur('syslin.Inverse', 'A must be an array of 2-axis')
    if n != m:
        ctrl.erreur('syslin.Inverse', 'Non square matrices cannot be inverted')
    if x.shape[0] != n or b.shape[0] != n:
        ctrl.erreur('syslin.Inverse', f'Non compatible input (dim(x)={x.shape}, dim(b)={x.shape}).')
    
    # Initializing Vars
    r = np.zeros(n)
    
    # Computing Residual
    for i in range(n):
        for j in range(n):
            r[i] += A[i,j]*x[j]
        r[i] -= b[i]
    # r /= n
    
    return r
# ============================================================================
def Residu_tridiag(L, D, U, b, x):
    n = len(D)
    r = np.zeros(n)
    r[0] = D[0]*x[0] + U[0]*x[1] - b[0]
    r[n-1] = L[n-1]*x[n-2] + D[n-1]*x[n-1] - b[n-1] 
    for i in range(1,n-1):
        r[i] = L[i]*x[i-1] + D[i]*x[i] + U[i]*x[i+1] - b[i]
    
    return r
# ============================================================================
def Norme2(r):
    """Compute the second norm of a vector r
    
    >>> Norme2([3,4]) # sqrt(3*3 + 4*4)
    >>> 5.0
    """
    norm = 0.
    
    for i in range(len(r)):
        norm += r[i]*r[i]
    norm = np.sqrt(norm)
    
    return norm
# ============================================================================
def Pivotation(A, n, p, b=None):
    """
    Find the maximum value in column p and switch rows in order to move it to position (p,p).
    
    **This method modifies A and b directly**

    Parameters
    ----------
    A : ndarray
         Coefficient Matrix.
    b : ndarray (optional)
        Right hand side of the linear system.
    n : int 
        System's order (number of unkowns).
    p : int
        Step number of Gauss-Jordan Elimination.
    
    Returns
    -------
        True if a row swapping was done, else None.
    
    """
    # Look for max value in column p
    k = p + np.argmax(abs(A[p:,p]))
    
    # Swapping rows
    if k != p: # using the power of numpy ndarrays, useful for larger arrays
        A[[p,k]] = A[[k,p]]
        if b is not None: b[[p,k]] = b[[k,p]]
        return True
# ============================================================================
def Inverse2(A):
    """
    Compute the (multiplicative) inverse matrix of A.

    Parameters
    ----------
    A : array_like
        Martix to be inversed.

    Returns
    -------
    invA : ndarray
        The inverse matrix of A.

    """
    # Input Control
    A = check_ndarray(A)
    try:
        n, m = A.shape
    except ValueError:
        ctrl.erreur('syslin.Inverse2', 'A must be an array of 2-axis')
    if n != m:
        ctrl.erreur('syslin.Inverse2', 'Non square matrices cannot be inverted')
    
    # Initializing Vars
    invA = np.zeros_like(A)
    c = np.zeros(n)
    c[n-1] = 1.
    
    # Computing 
    for i in range(n):
        c[i],c[i-1] = c[i-1],c[i]
        c = Gauss_Jordan(A, c)
        invA[i] = c
        
    return invA.T
# ============================================================================
def check_ndarray(*args, silent=True):
    """
    Check if given arguments are of type 'numpy.ndarray', else convert them.

    Parameters
    ----------
    *args : any_type
        The objects to be verified.
    silent : boolean, optional
        Do not to display the message 'converting array'.
        
    Raises
    ------
    TypeError
        One of the obejcts are not array_like.
    
    Returns
    -------
    A : numpy.ndarray
        The ndarray object.
    """
    ret = []
    for A in args:
        if not isinstance(A, np.ndarray):
            try:
                if not silent: print(f"Converting argument {args.index(A)} to ndarray...")
                A = np.array(A, dtype=np.float64)
            except:
                raise TypeError(f"Argument in position {args.index(A)} cannot be converted to ndarray")
        ret.append(A)
    return ret[0] if (len(ret) <= 1) else tuple(ret)
# ==============================================================