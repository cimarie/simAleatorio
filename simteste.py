
# coding: utf-8

# In[7]:

### Testes para a simulacao ###

# Valores padrao: N=5000, n=26, pA=1., b=0., c=10., deltaf=0.01, deltac=0.01, mutacao=0.0001, alpha=2., beta=0., pmig=0. #

from joblib import Parallel, delayed
from retest import simula, initSim
import multiprocessing as mp
import numpy as np

# N=1, n=10000, pA=1., b=0., c=10., deltaf=0.01, deltac=0.01, mutacao=0.0001, alpha=2., beta=0., pmig=0. #
N,n,pA,b,c,deltaf,mutacao,alpha,beta,pmig=1,10000,1.,0.,10.,0.01,0.0001,2.,0.,0.
A = pA*N*n
grupos,lfit,lfit_m,mpvencer = initSim(N,n,A,b,c,deltaf,alpha)

def bobo(i):
    return simula(N,n,mutacao,beta,pmig,grupos,lfit,lfit_m,mpvencer,i)

nc = 8 
l_its = Parallel(n_jobs=nc)(delayed(bobo)(i) for i in range(20))

# Mostra quantas simulacoes acabaram antes do maximo e a media de iteracoes que elas levaram
v = np.array(l_its)<5000
print "De 20 simulacoes, %d acabaram antes do maximo IT=5000" %(np.count_nonzero(v))
if np.count_nonzero(v)>0:
    w = np.take(l_its, np.nonzero(v)[0])
    media = sum(w)/len(w)
    print "Em media, a populacao de altruistas acaba em %.3f iteracoes" %media

