# -*- coding: utf-8 -*-
import numpy as np
import scipy.misc
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import pylab
import eigenvalue as eg
import numpy.linalg as LA
from math import *
from fitness import *
from selmig import *
from c_critico import *

def imprime_figura(data1, data2, label, tipo, titulo, nome):

    if tipo == 'av' or tipo == 'mc':
        plt.plot(data1,data2[0,:],label=label[0])
        for i in range(len(label)-1):
            plt.hold(True)
            plt.plot(data1,data2[i+1,:],label=label[i+1])
        plt.title(titulo)
        
        if tipo == 'av':
            plt.xlabel("m")
            plt.ylabel("Autovalor")

            pylab.xlim([0.0,0.6])
            pylab.xticks(np.arange(0.0,0.6,0.1))
            pylab.ylim([0.0,1.8])
            pylab.yticks(np.arange(0.0,1.8,0.2))

        elif tipo == 'mc':
            w, h = plt.figaspect(1.2)
            plt.figure(figsize=(w,h),dpi=100)
            plt.xlabel("delta")
            plt.ylabel("m critico")
            plt.autoscale(True)

        pylab.legend(loc='upper right')
        plt.hold(False)

    elif tipo == 'bm':
        w, h = plt.figaspect(1)
        plt.figure(figsize=(w,h),dpi=100)

        plt.plot(data1, data2)

        plt.xlabel("beta")
        plt.ylabel("m_critico")
        plt.title(titulo)
        
        pylab.xlim([-0.1,1.1])
        pylab.xticks(np.arange(0.0,1.05,0.2))
        pylab.ylim([-0.1,1.1])
        pylab.yticks(np.arange(0.0,1.05,0.2))

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(nome+".png")

    plt.close()


def salva_txt(data1, data2, tipo, nome):
    arq = open(nome+".txt",'w')
    if tipo == 'av' or tipo == 'mc':
        for cont in range(len(data1)):
            arq.write(str(data1[cont]))
            for i in range(data2.shape[0]):
                arq.write("\t")
                arq.write(str(data2[i,cont]))
            arq.write("\n")

    elif tipo == 'bm':
        for cont in range(len(data2)):
            arq.write(str(data1[cont]))
            arq.write("\t")
            arq.write(str(data2[cont]))
            arq.write("\n")

    arq.close()

# Cria figura autovalores vs. m para diferentes deltas
def autovalores_m(w0,b,c,n,beta,opcao):

    # Cria vetor com os valores de migracao que serao usados
    vetor_mig = np.arange(0,1.01,0.01,dtype=float)

    # Abre arquivo com os deltas
    with open('dados.txt', 'r') as f:
        exemplos = [float(x) for x in f.readline().split()]

    # Cria vetor para receber os autovalores da matriz
    autovalores = np.empty((len(exemplos),len(vetor_mig)))

    # Cria vetor para receber labels
    label = []

    for i in range(len(vetor_mig)):
        for j in range(len(exemplos)):
            m = vetor_mig[i]
            deltaf = deltac = exemplos[j]
            # label para colocar no grafico
            label.append('delta = ' + str(delta) + ' m_c = %.3f' %(m_critico(w0,b,c,n,deltaf,deltac,beta)))

            S,M = initSM(n,m,deltaf,deltac,b,c,w0,beta)
            autovalores[j][i] = eg.av_dominante(np.dot(M,S))

    #plt.figure(figsize=(30.6,3.48), dpi=100)
    
    nome = "bowles_delta_beta="+str(beta)
    if opcao==0:
        titulo = "Autovalor ("+"b = "+str(b)+", c = "+str(c) + ", n = " + str(n) + ", beta = " + str(beta) + ")"
        imprime_figura(vetor_mig,autovalores,label,titulo,'av',nome)

    else:
        salva_txt(vetor_mig,autovalores,'av',nome)

# Funcao rho: calcula (autovalor dominante - 1) a partir dos parametros
def rho(m, w0, b, c, n, deltaf, deltac , beta):
    S,M = initSM(n,m,deltaf,deltac,b,c,w0,beta)
    return eg.av_dominante(np.dot(M,S))-1

# Encontra m critico ate o qual todos os autovalores sao maiores que 1 (se 0<m<1, usa o metodo da secante para achar raizes)
def m_critico(w0,b,c,n,deltaf, deltac,beta):

    fa = rho(0.,w0,b,c,n,deltaf, deltac,beta)
    fb = rho(1.,w0,b,c,n,deltaf, deltac ,beta)

    if(fa*fb>0):
        if(fa>0):
            return 1.
        else:
            return 0.
    # brenth encontra o zero da funcao rho no intervalo [0.,1.]
    return brentq(rho, 0., 1., xtol=0.001, args=(w0,b,c,n,deltaf, deltac,beta))

# Migracao critica pelo metodo de interpolacao
def m_critico2(w0,b,c,n,deltaf,deltac,beta):
    precisao = 0.00001

    fa = rho(0.,w0,b,c,n,deltaf, deltac,beta)
    fb = rho(1.,w0,b,c,n,deltaf, deltac ,beta)

    if(fa*fb>0):
        if(fa>0):
            return 1.
        else:
            return 0.

    vec_m = [0.,1.]
    fp = map(lambda i: rho(i,w0,b,c,n,deltaf,deltac,beta), vec_m)
    npt = 100
    ind = 0 
    while ind < len(fp)-1 and fp[ind]*fp[ind+1] < 0:
        vec_m = np.arange(vec_m[ind],vec_m[ind+1]+1./npt,1./npt)
        fp = map(lambda i: rho(i,w0,b,c,n,deltaf,deltac,beta), vec_m)
        fp = np.array(fp)
        ind = np.where(fp>0)[0][-1]
        npt = 10*npt

    return vec_m[ind] 

def encontra_mc(w0,beta,b,c,opcao):

    # Abre arquivo com os tamanhos dos grupos
    with open('ndados', 'r') as f:
        exemplos = [int(x) for x in f.readline().split()]

    vetor_delta = np.arange(0.01,1.0,0.01,dtype=float)
    vetor_migc = np.empty((len(exemplos),len(vetor_delta)))
    label = []

    for i in range(len(exemplos)):
        n = exemplos[i]
        print n
        label.append('n = ' + str(n))

        for j in range(len(vetor_delta)):
            deltac = deltaf = vetor_delta[j]
            vetor_migc[i][j] = m_critico(w0,b,c,n,deltaf,deltac,beta)
   
    nome = "bowles_m_critico_beta="+str(beta)

    if opcao == 0:
        titulo = "Migracao critica ("+"b = "+str(b)+", c = "+str(c) + ", beta = " + str(beta) + ")"
        imprime_figura(vetor_delta,vetor_migc,label,'mc',titulo,nome)

    else:
        salva_txt(vetor_delta,vetor_migc,'mc',nome)

# Gera figura beta(m) sem as cores ciano magenta
def beta_m(w0, b, c, n, deltaf, deltac ,opcao):

    r = 5000
    vetor_beta = np.arange(0.0,1.+1./r,1./r,dtype=float)
    vetor_m = np.empty(len(vetor_beta))

    for i in range(len(vetor_beta)):
        vetor_m[i]=m_critico(w0,b,c,n,deltaf, deltac ,vetor_beta[i])
        
    nome = "mcritico_versus_beta_delta="+str(deltaf)
    if opcao==0:

        titulo = "m_critico x beta ("+"b = "+str(b)+", c = "+str(c) + ", n = " + str(n) + ",\ndeltaf = " + str(deltaf) + ", deltac = " + str(deltac) + ")"
        imprime_figura(vetor_beta,vetor_m,[],'bm',titulo,nome)

    else:
        salva_txt(vetor_beta,vetor_m,'bm',nome)

# gera figura beta(m) divisao nas cores ciano e magenta
def beta_m2(w0, b, c, n, deltaf, deltac):

    r = 100

    vetor_beta = []
    vetor_m = []

    for i in range(r+1):

        vetor_beta = np.append(vetor_beta,np.arange(0.0,1.+1./r,1./r,dtype=float))

        vetor_m = np.append(vetor_m,[float(i)/r]*(r+1))

    print len(vetor_beta)

    for i in range(len(vetor_beta)):

        av = rho(vetor_m[i],w0,b,c,n,deltaf,deltac,vetor_beta[i])+1

        z[i] = 'c' if av<1 else 'm'

    w, h = plt.figaspect(1)
    plt.figure(figsize=(w,h),dpi=100)

    # Mostra um grafico a migracao critica vs. delta
    plt.scatter(vetor_m,vetor_beta,c=z,edgecolors='none')

    # Label
    plt.xlabel("m")
    plt.ylabel("beta")

    # Titulo
    titulo = "beta x m ("+"b = "+str(b)+", c = "+str(c) + ", n = " + str(n) + ", delta = " + str(deltaf) + ")"
    plt.title(titulo)

    # Limite e intervalos do eixo x
    pylab.xlim([0.0,1.0])
    pylab.xticks(np.arange(0.0,1.1,0.1))

    # Limite e intervalos do eixo y
    pylab.ylim([0.0,1.0])
    pylab.yticks(np.arange(0.0,1.1,0.1))

    plt.tight_layout()

    # Salva fig
    nome = "bowles_m_beta_delta="+str(deltaf) 
    plt.savefig(nome+".png")

    # Fecha fig
    plt.close()


# Faz o gráfico em barras de autovetores
def autovetores(n,b,c,w0,beta):

    vetor_mig = np.arange(0.02,0.06,0.01,dtype=float)
    vetor_delta = np.array((0.0001,0.001,0.01),dtype=float)

    avec = np.empty((len(vetor_mig),n))

    j = 0

    vetor = np.arange(0.8,n,1,dtype=float)

    # Define as caracteristicas da figura gerada
    plt.figure(1)

    plt.rc('text',usetex=True)
    plt.rc('font', family='serif')

    #plt.tight_layout()

    f, ax = plt.subplots(3, 4,sharey=True)

    aux = np.arange(n)+1

    for delta in vetor_delta:
        deltaf = delta
        deltac = delta

        # Cria matriz de selecao
        S = matriz_selecao(n,deltaf,deltac,b,c,w0,beta)
        i = 0

        for m in vetor_mig:

            # Cria matriz de migracao
            M = matriz_migracao(n,m)

            # Calcula autovetores
            avec[i] = eg.calcula_evec(np.dot(S,M))

            avec[i] = aux*avec[i]/(np.dot(aux,avec[i]))

            #avec[i] = avec[i]/sum(avec[i])

            #ax=plt.subplot(3,4,4*j+i)

            ax[j,i].set_xticks((1,5,10,15,20))

            if i==0:
                ax[j,0].set_ylabel(r'$\nu$')

            ax[j,i].set_yticks(np.arange(0.0,0.5,0.1))

            ax[j,i].autoscale_view(True,True,True)

            ax[j,i].set_xlabel(r"$k$")

            lab = r'$\rho$='+"%.3f" % (eg.av_dominante(np.dot(S,M)))+"\n"+r'$m_c$='+"%.3f" % (m_critico(w0,b,c,n,deltaf,deltac,beta))+", "+r'$\delta$='+str(deltaf)

            ax[j,i].set_title(lab,fontsize=8)

            plt.grid()

            # Atualiza contador
            i = i+1

            print 4*j+i

            ax[j,i-1].bar(vetor,avec[i-1]/sum(avec[i-1]),0.4)

            plt.legend(loc="center")

        j = j+1

    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.3, hspace=0.7)

    # titulo da figura
    titulo = "Equilibrium distribution (b="+str(b)+", c="+str(c)+", n="+str(n)+","+r'$\beta$='+str(beta)+")"
    plt.suptitle(titulo) 

    # Salva fig
    nome = "b_eigenvec_b="+str(b)+"c=1.eps"
    plt.savefig(nome)

    # Fecha fig
    plt.close()

## Main
def main():

    w0 = 1.             # fitness basal
    beta = 0.           # probabilidade de combate
    b = 0.              # beneficio gerado por altruista
    c = 10              # custo para um altruista
    n = 26
    delta = 0.01

    print m_critico2(w0,6,3,n,delta,delta,beta)
#
##    #m = 0.1
##
##    deltaf = delta
##    deltac = delta
##
##    # Cria matriz de selecao
##    #S = matriz_selecao(n,deltaf,deltac,b,c,w0,beta)
#
#    # Cria matriz de migracao
##    M = matriz_migracao(n,m)
#
##    vetor = eg.calcula_evec(np.dot(S,M))
#
#    # Para Fst = 0.04 => m = 0.8
#    m = 0.8
#    deltac = .01
#    deltaf = delta
#    #b=5*c
#    
#   # for beta in [0.414, 0.2, 0.09]:
#   #     encontra_cc(w0,beta,b,m,0)
#
#    #for deltaf in [0.01,0.5,0.8]:
#    #    deltac = deltaf
#    #    beta_c(w0, b, m, n, deltaf,deltac,0)
#    #    beta_m(w0, b, c, n, deltaf,deltac,0)
#
#    #deltaf = 0.01
#    #deltac = 0.01927
#
#    #for beta in [0.414, 0.2, 0.09]:
#    #    print beta/2
#    #    print c_critico(w0,b,m,n,deltaf,deltac,beta)
#
#    #autovalores_m(w0,b,c,n,beta,1)
#    #autovalores_m(w0,b,c,n,beta,0)
#
#    #beta_c(w0,b,m,n,deltaf,deltac,0)
#
#    #autovetores(n,b,c,w0,beta)
#
#    encontra_mc(w0,beta,b,c,0)
#
#
#
if  __name__ =='__main__':main()
