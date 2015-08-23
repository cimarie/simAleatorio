#coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import pylab
from c_critico import *

pylab.rcParams['text.usetex']=True
pylab.rcParams['font.family'] = 'serif'
pylab.rcParams['font.serif'] = 'cm'
pylab.rcParams['font.size'] = 16

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
    resbow = []
    resmeu = []
    vbeta = []
    for mort in [0.207, 0.100, 0.045]:
        resbow.append(devolve_cmax(fst,mort))
        beta = 2*mort
        vbeta.append(beta)
        resmeu.append(c_critico(w0,b,m,n,delta,alpha,beta))

    params = {"b":b, "n": n, "m": m, "\\alpha": alpha, "\\delta": delta}
    str_params = "\n"
    counter = 0
    lista = list(params.keys())
    lista.sort()
    for nome_param in lista:
        counter += 1
        str_params += r"$%s = $%s, " %(nome_param, str(params[nome_param]))
        if counter > 4 and nome_param != params.keys()[-1]:
            str_params += "\n"
    titulo = ""
    titulo = titulo + str_params[:-2]

    w, h = plt.figaspect(1)
    plt.figure(figsize=(w,h), dpi=300)

    lw = 2.0
    plt.plot(data1, data2, linewidth=lw)

    plt.grid(True)
    plt.tight_layout()
    plt.autoscale(True)

    plt.title(titulo)
    plt.savefig("max_cost_delta=%s.tif")
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
