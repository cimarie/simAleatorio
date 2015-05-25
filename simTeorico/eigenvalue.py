# -*- coding: utf-8 -*-
#eigenvalue.py

import scipy.linalg
import numpy as np

# Calcula autovalor dominante
def av_dominante(Matriz):

    Matriz = np.delete(Matriz,0,0)

    Matriz = np.delete(Matriz,0,1)

    av = scipy.linalg.eigvals(Matriz)

    x = np.argsort(-abs(av))

    return np.real(av[x[0]])

# Calcula vetor associado ao autovalor dominante
def calcula_evec(Matriz):

    n = Matriz.shape[0]-1

    Matriz = np.delete(Matriz,0,0)

    Matriz = np.delete(Matriz,0,1)

    av,avec = scipy.linalg.eig(Matriz,left=True,right=False)

    x = np.argsort(-abs(av))

    return np.fabs(avec[:,x[0]])

# Devolve a soma da coluna i
def soma_col(M,i):

    return sum(M[:,i])

# Devolve vetor com o valor maximo e minimo, considerando as somas de cada coluna (para estimativa de autovalor dominante)
def maxmin_cols(M,n):

    aux = np.arange(n+1)

    v = soma_col(M,aux)

    v = np.sort(v)

    return np.array([v[0],v[n]])
