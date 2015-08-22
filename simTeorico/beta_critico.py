from selmig import *
import eigenvalue as eg
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import pylab

# Funcao rhob: calcula (autovalor dominante - 1) a partir dos parametros
def rhob(beta,w0,b,c,m,n,delta,alpha):
    
    S,M = initSM(n, m, delta, alpha, b, c, w0, beta)
    return eg.av_dominante(np.dot(M,S))-1

# Encontra beta critico ate o qual todos os autovalores sao maiores que 1 (se 0<m<1, usa o metodo da secante para achar raizes)
def beta_critico(w0,b,c,m,n,delta,alpha):

    fa = rhob(0.,w0,b,c,m,n,delta,alpha)

    fb = rhob(1.,w0,b,c,m,n,delta,alpha)

    if(fa*fb>0):
        if(fa>0):
            return 1.
        else:
            return 0.

    # brenth encontra o zero da funcao rhob no intervalo 
    return brentq(rhob, 0., 1., xtol=0.001, args=(w0,b,c,m,n,delta,alpha))

# Faz curva beta x alpha com os valores de beta que suportam c = 0.03
def encontra_beta(w0, b, c, n, delta, m):

    f = open('alpha_beta_m=%.3f.txt'%m,'w')
    list_alpha = np.arange(0., 3.005, 0.01) 
    list_beta = []
    for alpha in list_alpha:
        beta = beta_critico(w0,b,c,m,n,delta,alpha)   
        list_beta.append(beta)
        f.write(str(beta))
        f.write("\t")
        f.write(str(alpha))
        f.write("\n")
    f.close()

# Funcao rhoal: calcula (autovalor dominante - 1) a partir dos parametros
def rhoal(alpha,w0,b,c,m,n,delta,beta):
    
    S,M = initSM(n, m, delta, alpha, b, c, w0, beta)
    return eg.av_dominante(np.dot(M,S))-1

# Encontra beta critico ate o qual todos os autovalores sao maiores que 1 (se 0<m<1, usa o metodo da secante para achar raizes)
def alpha_critico(w0,b,c,m,n,delta,beta):

    fa = rhoal(0.,w0,b,c,m,n,delta, beta)

    fb = rhoal(3.,w0,b,c,m,n,delta, beta)

    if(fa*fb>0):
        if(fa>0):
            return 1.
        else:
            return 0.

    # brenth encontra o zero da funcao rhoal no intervalo 
    return brentq(rhoal, 0., 3., xtol=0.001, args=(w0,b,c,m,n,delta,beta))

# Faz curva beta x alpha com os valores de beta que suportam c = 0.03
def encontra_alpha(w0, b, c, n, delta, m):

    f = open('beta_alpha_m=%.3f.txt'%m,'w')
    list_beta = np.arange(0., 0.5004, 0.01) 
    list_alpha = []
    for beta in list_beta:
        alpha = alpha_critico(w0,b,c,m,n,delta,beta)   
        f.write(str(beta))
        f.write("\t")
        f.write(str(alpha))
        f.write("\n")
        list_alpha.append(alpha)
    f.close()
