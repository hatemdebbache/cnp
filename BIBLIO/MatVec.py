"""
        Module MatVec
        
        Methods pour des operations matricielles et vectorielles
"""
import numpy as np
from math import sqrt

# ==============================================================
def vec_interval(a, b, npts):
    """ 
        Return an array of n equidistant points in interval [a,b]
    """
    x = np.zeros(npts)
    pas = (b-a)/(npts - 1)
    for i in range(npts):
        x[i] = a + i*pas
    return x
# ==============================================================
def modif_mat(x, v, i1, j1, i2, j2):
    """ 
        Modify all values in square of diagonal {(i1,j1), (i2,j2)} to be v
    """
    for i in range(i1,i2+1):
        for j in range(j1,j2+1):
            x[i,j] = v
# ==============================================================
def vec_valeur(n, v):
    """ 
        Return a vector filled with value v
    """
    vec = np.zeros(n)
    for i in range(n):
        vec[i] = v

    return vec
# ==============================================================
def vec_indice(n):
    """ 
        Return a vector where each component is equal to its index
        (equivalent to np.arange(n))
    """
    vec = np.zeros(n)
    for i in range(n):
        vec[i] = i 
        
    return vec
# ==============================================================
def inc_vec(x, m):
    """ 
        Add m to each component to vector x
    """
    for i in range(len(x)):
        x[i] += m
# ==============================================================
def copie_vec(x):
    """ 
        Return a copy of the input
    """
    y = np.zeros(x.shape)
    y += x
    
    return y
# ==============================================================
def module(u):
    """ 
        Return the 2-norm of a vector u (length of vector)
    """
    s = 0.
    for i in range(len(u)):
        s += u[i]*u[i]
        
    return sqrt(s)
# ==============================================================
def dot_v(u,v):
    """
        Return the scalar product of vectors u and v (u.v)
    """
    if len(u) != len(v):
        print("Vectors must have same dimension")
        return
    s = 0.
    for i in range(len(u)):
        s += u[i]*v[i]
    
    return s
# ==============================================================
def cross_v(u, v):
    """ 
        Returns the cross product of vectors u and v (u^v)
    """
    if len(u) != 3 or len(v) != 3:
        # raise DimensionError("Product doesn't exist for these vectors")
        pass
    w = np.zeros(3)
    
    w[0] = u[1]*v[2] - u[2]*v[1]
    w[1] = u[2]*v[0] - u[0]*v[0]
    w[2] = u[0]*v[1] - u[1]*v[0]
    
    return w 
# ==============================================================
def mat_triangulaire(n, m, v, sup = True):
    M = np.zeros((n,m))
    for i in range(n):
        if not sup:
            for j in range(i+1):
                M[i,j] = v
        else:
            for j in range(i,m):
                M[i,j] = v
    return M
# ==============================================================
def dot_m(A, B):
    """ 
        Return P the matrix product of A and B (must be of 2 Axis)
    """
    d1, d2 = A.shape, B.shape
    if len(d1) != 2 or len(d2) != 2:
        print("Please provide matrices with 2 axis")
    elif d1[0] != d2[1]:
        print("Matrix product not possible")
    else:
        n, m = d1[0], d2[1]
        ret = np.zeros((n, m))
        for i in range(n):
            for j in range(m):
                for p in range(d1[1]):
                    ret[i,j] += A[i,p]*B[p,j]
        
        return ret
# ==============================================================
def trans_m(M):
    """ 
        Return the transpose matrix of M
    """
    dim = list(M.shape)
    if len(dim) != 2:
        print("please provide a matrix with 2 axis")
        return
    m, n = dim
    ret = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            ret[i,j] = M[j,i]
            
    return ret
# ==============================================================
def copie_mat(A, i1, j1, i2, j2):
    """ 
        Return a copy of block [(i1,j1) ... (i2,j2)] of matrix A
    """
    di = i2 - i1 + 1
    dj = j2 - j1 + 1
    ret = np.zeros((di, dj))
    for i in range(di):
        for j in range(dj):
            ret[i,j] = A[i+i1,j+j1]
    return ret
# ==============================================================
def trace(A):
    """ 
        Return Trace(A) = Sum(Diag(A))
    """
    s = 0.
    dim = A.shape
    if len(dim) != 2:
        print("please provide a matrix with 2 axis")
        return
    n = min(dim)
    for i in range(n):
        s += A[i,i]
    
    return s
# ==============================================================
def is_magic(A):
    """ 
        Return whether A is a magic square (True) or not (False)
        
        Magic squre: n x n matrix of the distinct elements from 1 to n² where
            the sum of any row, column, or diagonal is always equal to the same number
        Ex: A = [[2 7 6]
                 [9 5 1]
                 [4 3 8]]
    """
    dim = A.shape
    if len(dim) != 2:
        print("please provide a matrix with 2 axis")
        return
    t = trace(A)
    A1 = trans_m(A)
    if t == trace(A1):
        def check(A, t):
            for i in range(dim[0]):
                s = somme_ligne(A, i)
                if t != s:
                    return False
            return True
        if check(A, t):
            return check(A1, t)
    return False
# ==============================================================
def somme_ligne(A, ligne):
    """ 
        Sum the elements of A[ligne]
    """
    s = 0.
    n = A.shape[1]
    for i in range(n):
        s += A[ligne, i]
    return s
# ==============================================================
def occurence(A, v):
    """ 
        Return the number of occurences of a value v in A
    """
    dim = A.shape
    if len(dim) != 2:
        print("please provide a matrix with 2 axis")
        return
    n, m = dim
    ret = 0
    for i in range(n):
        for j in range(m):
            if A[i,j] == v:
                ret += 1
    return ret
# ==============================================================
def create_mat(n, m):
    """ 
        Return a matix with each element equal to the sum of its ligne and column
    """
    x = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            x[i,j] = i+j
            
    return x
# ==============================================================