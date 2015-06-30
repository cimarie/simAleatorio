# -*- coding: utf-8 -*-
import numpy as np
import scipy.misc
import eigenvalue as eg
import numpy.linalg as LA
from itertools import chain, izip
from fitness import *

# Probabilidade de vencer - recebe k e n
# calcula a probabilidade de vencer um combate
def pVencer(alpha,k,n):
#    arg1 = -4*alpha*k/n
#    return 1/(1+np.exp(arg1))
    return 0.5 + alpha*k/n

# Select - recebe indice k,l e o numero total n
# calcula o numero medio de grupos do tipo k que viram tipo l por força de selecao
def sel(k,l,n,delta, alpha,b,c,w0,beta):
    w = fitness(0,k,n,delta,b,c,w0)
    wm = fitness_m(k,n,delta,b,c,w0)
    p = 0 if wm == 0 else k*w/(n*wm)
    termoInd = scipy.misc.comb(n,l)*p**l*(1-p)**(n-l)
    #termoGr = (1-beta+2*beta*pVencer(alpha,k,n)) if beta > 0 else wm
    termoGr = 1-beta+2*beta*pVencer(alpha,k,n)
    return termoGr*termoInd

# Matriz de selecao - recebe n e parametros delta,b,c,w0 e beta
def matriz_selecao(n,delta, alpha,b,c,w0,beta):
    #S = np.empty([n+1,n+1], dtype=float)
    grid = np.indices((n+1,n+1))
    #S = map(lambda i,j: sel(i,j,n,delta,alpha,b,c,w0,beta), grid[0], grid[1])
    indices = izip(chain.from_iterable(grid[0]),chain.from_iterable(grid[1]))
    S = np.array([sel(k,l,n,delta,alpha,b,c,w0,beta) for k,l in indices])
    S = S.reshape(n+1,n+1)
    #for elem in indices:
    #    k, l = elem
    #    S[k][l] = sel(k,l,n,delta,alpha,b,c,w0,beta)

    return S

# Mig - recebe indice k,l e o numero total n
# Calcula probabilidade de um grupo k virar tipo l em relacao ao efeito de migracao
def mig(k,l,m):

    return scipy.misc.comb(k,l)*((1-m)**l)*(m**(k-l))

# Matriz do numero medio de grupos que eram do tipo k e viram tipo l depois da migracao
def matriz_mig1(n,m):

    M = np.fromfunction(lambda i,j: mig(i,j,m),(n+1,n+1))

    return M

# Calcula o numero medio de grupos que recebem um altruista (eram do tipo 0 e viram do tipo 1)
def func(k,m):

    return m*k

# Matriz do numero medio de grupos que recebem altruistas
def matriz_mig2(n,m):

    T = np.zeros((n+1,n+1),dtype=float)

    T[:,1] = func(np.arange(n+1),m)

    return T

# Matriz de migracao - recebe n
def matriz_migracao(n,m):


    M = np.identity(n+1) if m==0 else matriz_mig1(n,m)

    T = matriz_mig2(n,m)

    return M+T

# Inicializa matrizes de selecao
def initSM(n, m, deltaf, alpha, b, c, w0, beta):
    
    # Cria matriz de selecao
    S = matriz_selecao(n,deltaf,alpha,b,c,w0,beta)

    # Cria matriz de migracao
    M = matriz_migracao(n,m)

    return S, M
