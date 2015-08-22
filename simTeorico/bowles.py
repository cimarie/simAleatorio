#coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import pylab
from c_critico import *

def devolve_lambda(mu, Fst, b, c, n):

    return (c-b/n)/(4*mu*((Fst/(1-Fst))+(1./n)))

def devolve_cmax(fst, mortalidade):
    n = 26
    lambdaA = 2.
    return 4*lambdaA*mortalidade*((fst/(1-fst))+1./n)

def calcula_r(n, fst):

    return (n*fst-1)/(n-1)

def calcula_m(n, r):
    A = 1/(1+n*r/(1-r)) 
    print 1-A    
    if 1-A > 0:
        m =  1+sqrt(1-A), 1-sqrt(1-A)
    else:
        m =  0, 0
    if max(m) > 1:
        if min(m) > 0:
            return min(m)
        else:
            return 0.
    else: 
        return max(m)

def main():
    
    n = 26
    fst = 0.040
    r = calcula_r(n,fst)
    print u"\nO valor de r é: %.7f\n" %r
    m = calcula_m(n,r)
    print u"\nO valor de m é:%.7f" %m
    w0, b, delta, alpha = 1, 0., 0.1, 2.
    for mort in [0.207, 0.100, 0.045]:
        print devolve_cmax(fst,mort)
        beta = 2*mort
        print c_critico(w0,b,m,n,delta,alpha,beta)


    #c = 3
    #b = 5 
    #Fst_vec = [0.022,0.075,0.170]
    #n = 26 
    #label = []
    #mu_vec = np.arange(0.0001,0.26,0.00001, dtype=float)
    #val_lambda = np.empty((len(Fst_vec), len(mu_vec)))

    #for i in range(len(Fst_vec)):
   
    #    label.append("Fst = " + str(Fst_vec[i]))
    #
    #    for j in range(len(mu_vec)):
    #        val_lambda[i][j]=devolve_lambda(mu_vec[j],Fst_vec[i],b,c,n)

    #plt.plot(mu_vec,val_lambda[0,:],label=label[0])
    #plt.hold(True)
    #plt.plot(mu_vec,val_lambda[1,:],label=label[1])
    #plt.hold(True)
    #plt.plot(mu_vec,val_lambda[2,:],label=label[2])

    #plt.xlabel(r'$\delta$')
    #plt.ylabel(r'$\lambda_A$')

    #plt.grid(True)

    #pylab.legend()
    #plt.hold(False)

    #pylab.xlim([0.0,0.05])
    #pylab.xticks(np.arange(0.0,0.26,0.05))

    #pylab.ylim([0.0,3.0])
    #pylab.yticks(np.arange(0.0,3.1,0.5))

#    plt.autoscale(True)

    #plt.tight_layout()

    #nome = "figura_bowles_c.png"
    #plt.savefig(nome)

    #plt.close()

if __name__ == '__main__':
    main()
