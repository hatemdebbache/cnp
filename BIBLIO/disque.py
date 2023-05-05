# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:14:56 2020

@author: user
"""
import ctrl
import numpy as np
# ==============================================================
def LireVec2(nom_fichier):
    x, y = [], []
    with open(nom_fichier, 'r') as f:
        for line in f.readlines():
            s = line.split()
            x.append(float(s[0]))
            y.append(float(s[1]))
    npts = len(x)
    return npts, x, y
# ==============================================================
def LireTridiag(nom_fichier):

    L, D, U, b = [],[],[],[]
    with open(nom_fichier, 'r') as f:
        for line in f.readlines():
            s = line.split()
            L.append(float(s[0]))
            D.append(float(s[1]))
            U.append(float(s[2]))
            b.append(float(s[3]))
    npts = len(L)
    
    return npts, L, D, U, b
# ==============================================================
def LireMat(nom_fichier):

    with open(nom_fichier, 'r') as f:
        header = f.readline().split()
        N, M = int(header[0]), int(header[1])
        A = [[float(s) for s in line.split()] for line in f.readlines()]
        A = np.array(A)
        assert(A.shape == (N, M))

    return N, M, A
# ==============================================================
def LireSys(nom_fichier):
    with open(nom_fichier, 'r') as f:
        header = f.readline().split()
        N, M = int(header[0]), int(header[1])
        X = [[float(s) for s in line.split()] for line in f.readlines()]
        X = np.array(X)
        assert(X.shape == (N, M + 1))
        A = X[:,:-1]
        b = X[:,-1]
    
    return N, M, A, b
# ==============================================================
def EcrireMat(M, nom_fichier, cs): 
    """
    Ecrire une matrice M dans un fichier text.
    
    Parameters
    ----------
    M : ndarray
        La matrice à écrire.
    nom_fichier : str
        Le nom de fichier de sortie, le chemin du fichier peut être inclut.
        (N'oblier pas d'ajouter .txt à la fin du nom)
    cs : int
        Le nombre de chiffres significatifs.
    
    Notes
    -----
    Be careful with the truncation effect when assigning 'cs'
    """
    try:
        n,m = M.shape
    except:
        ctrl.erreur('EcrireMat', 'M n\'est pas une matrice (de type ndarray)')
    
    with open(nom_fichier, 'w') as sortie:
        sortie.write(str(n) + ' ' + str(m))
        for i in range(n):
            sortie.write('\n')
            for j in range(m):
                sortie.write(str(format(M[i,j], f'.{cs}g')) + ' ')
# ==============================================================