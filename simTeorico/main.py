# -*- coding: utf-8 -*-
from joblib import Parallel, delayed
import os.path
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
from beta_critico import *

pylab.rcParams['text.usetex']=True
pylab.rcParams['font.family'] = 'serif'
pylab.rcParams['font.serif'] = 'cm'

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
            delta = exemplos[j]
            # label para colocar no grafico
            label.append('delta = ' + str(delta) + ' m_c = %.3f' %(m_critico(w0,b,c,n,delta,alpha,beta)))

            S,M = initSM(n,m,delta,alpha,b,c,w0,beta)
            autovalores[j][i] = eg.av_dominante(np.dot(M,S))

    #plt.figure(figsize=(30.6,3.48), dpi=100)
    
    nome = "bowles_delta_beta="+str(beta)
    if opcao==0:
        titulo = "Autovalor ("+"b = "+str(b)+", c = "+str(c) + ", n = " + str(n) + ", beta = " + str(beta) + ")"
        imprime_figura(vetor_mig,autovalores,label,titulo,'av',nome)

    else:
        salva_txt(vetor_mig,autovalores,'av',nome)

# Funcao rho: calcula (autovalor dominante - 1) a partir dos parametros
def rho(m, w0, b, c, n, delta, alpha , beta):
    S,M = initSM(n,m,delta,alpha,b,c,w0,beta)
    return eg.av_dominante(np.dot(M,S))-1

# Encontra m critico ate o qual todos os autovalores sao maiores que 1 (se 0<m<1, usa o metodo da secante para achar raizes)
def m_critico(w0,b,c,n,delta, alpha,beta):

    fa = rho(0.,w0,b,c,n,delta, alpha,beta)
    fb = rho(1.,w0,b,c,n,delta, alpha ,beta)

    if(fa*fb>0):
        if(fa>0):
            return 1.
        else:
            return 0.
    # brenth encontra o zero da funcao rho no intervalo [0.,1.]
    return brentq(rho, 0., 1., xtol=0.001, args=(w0,b,c,n,delta,alpha,beta))

# Migracao critica pelo metodo de interpolacao
def m_critico2(w0,b,c,n,delta,alpha,beta):
    precisao = 0.00001

    fa = rho(0.,w0,b,c,n,delta,alpha,beta)
    fb = rho(1.,w0,b,c,n,delta,alpha ,beta)

    if(fa*fb>0):
        if(fa>0):
            return 1.
        else:
            return 0.

    vec_m = [0.,1.]
    fp = map(lambda i: rho(i,w0,b,c,n,delta,alpha,beta), vec_m)
    npt = 100
    ind = 0 
    while ind < len(fp)-1 and fp[ind]*fp[ind+1] < 0:
        vec_m = np.arange(vec_m[ind],vec_m[ind+1]+1./npt,1./npt)
        fp = map(lambda i: rho(i,w0,b,c,n,delta,alpha,beta), vec_m)
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
            delta = vetor_delta[j]
            vetor_migc[i][j] = m_critico(w0,b,c,n,delta,alpha,beta)
   
    nome = "bowles_m_critico_beta="+str(beta)

    if opcao == 0:
        titulo = "Migracao critica ("+"b = "+str(b)+", c = "+str(c) + ", beta = " + str(beta) + ")"
        imprime_figura(vetor_delta,vetor_migc,label,'mc',titulo,nome)

    else:
        salva_txt(vetor_delta,vetor_migc,'mc',nome)

# Gera figura beta(m) sem as cores ciano magenta
def beta_m(w0, b, c, n, delta, alpha ,opcao):

    r = 5000
    vetor_beta = np.arange(0.0,1.+1./r,1./r,dtype=float)
    vetor_m = np.empty(len(vetor_beta))

    for i in range(len(vetor_beta)):
        vetor_m[i]=m_critico(w0,b,c,n,delta,alpha,vetor_beta[i])
        
    nome = "mcritico_versus_beta_delta="+str(delta)
    if opcao==0:

        titulo = "m_critico x beta ("+"b = "+str(b)+", c = "+str(c) + ", n = " + str(n) + ",\ndelta = " + str(deltaf) + ", alpha = " + str(alpha) + ")"
        imprime_figura(vetor_beta,vetor_m,[],'bm',titulo,nome)

    else:
        salva_txt(vetor_beta,vetor_m,'bm',nome)

# gera figura beta(m) divisao nas cores ciano e magenta
def beta_m2(w0, b, c, n, delta, alpha):

    r = 100

    vetor_beta = []
    vetor_m = []

    for i in range(r+1):

        vetor_beta = np.append(vetor_beta,np.arange(0.0,1.+1./r,1./r,dtype=float))

        vetor_m = np.append(vetor_m,[float(i)/r]*(r+1))

    print len(vetor_beta)

    for i in range(len(vetor_beta)):

        av = rho(vetor_m[i],w0,b,c,n,delta,alpha,vetor_beta[i])+1

        z[i] = 'c' if av<1 else 'm'

    w, h = plt.figaspect(1)
    plt.figure(figsize=(w,h),dpi=100)

    # Mostra um grafico a migracao critica vs. delta
    plt.scatter(vetor_m,vetor_beta,c=z,edgecolors='none')

    # Label
    plt.xlabel("m")
    plt.ylabel("beta")

    # Titulo
    titulo = "beta x m ("+"b = "+str(b)+", c = "+str(c) + ", n = " + str(n) + ", delta = " + str(delta) + ")"
    plt.title(titulo)

    # Limite e intervalos do eixo x
    pylab.xlim([0.0,1.0])
    pylab.xticks(np.arange(0.0,1.1,0.1))

    # Limite e intervalos do eixo y
    pylab.ylim([0.0,1.0])
    pylab.yticks(np.arange(0.0,1.1,0.1))

    plt.tight_layout()

    # Salva fig
    nome = "bowles_m_beta_delta="+str(delta) 
    plt.savefig(nome+".png")

    # Fecha fig
    plt.close()


# Faz o gráfico em barras de autovetores
def autovetores(n,b,c,w0,alpha,beta):

    #vetor_mig = np.arange(0.02,0.06,0.01,dtype=float)
    vetor_mig = [0.01, 0.1, 0.2, 0.5]
    vetor_delta = np.array((0.0001,0.001,0.01),dtype=float)

    avec = np.empty((len(vetor_mig),n))

    j = 0

    vetor = np.arange(0.8,n,1,dtype=float)

    # Define as caracteristicas da figura gerada
    #plt.figure(figsize=(6,6), dpi=300)

    #plt.rc('text',usetex=True)
    #plt.rc('font', family='serif')

    fig, ax = plt.subplots(nrows=3,ncols=4, sharey=True, dpi=300)
    #f, ax = plt.subplots(3, 4,sharey=True)
    plt.setp(ax, xticks=[1,5,10,15,20,25], yticks=np.arange(0.0,0.7,0.1))

    #fig.set_canvas(plt.gcf().canvas)

    aux = np.arange(n)+1

    for delta in vetor_delta:

        # Cria matriz de selecao
        S = matriz_selecao(n,delta,alpha,b,c,w0,beta)
        i = 0

        for m in vetor_mig:

            # Cria matriz de migracao
            M = matriz_migracao(n,m)

            # Calcula autovetores
            avec = eg.calcula_evec(np.dot(S,M))

            avec = aux*avec/(np.dot(aux,avec))

            avec = avec/sum(avec)

            print i, j
            print avec

            #ax=plt.subplot(3,4,4*j+i)

            #ax[j,i].set_xticks((1,5,10,15,20))

            if i==0:
                ax[j,0].set_ylabel(r'$\nu$', fontsize=8)

            #ax[j,i].set_yticks(np.arange(0.0,0.5,0.1))

            ax[j,i].autoscale_view(True,True,True)

            ax[j,i].set_xlabel(r"$k$", fontsize=8)

            texto = r'$\rho$='+"%.3f" % (eg.av_dominante(np.dot(S,M)))+"\n"+r'$m_c$='+"%.3f" % (m_critico(w0,b,c,n,delta,alpha,beta))

            plt.text(0.5,0.5,texto,horizontalalignment='center',\
                    verticalalignment='center',
                    transform=ax[j,i].transAxes,
                    bbox=dict(facecolor='white', alpha=0.5),
                    fontsize=8)
            lab = r'$\delta=$'+str(delta)+', '+r'$m=$'+str(m)
            ax[j,i].set_title(lab,fontsize=9)

            plt.grid()

            # Atualiza contador
            i = i+1

            ax[j,i-1].bar(vetor,avec,width=0.4)

            #plt.legend(loc="center")

        j = j+1

    plt.subplots_adjust(left=0.125, bottom=0.001, right=0.9, top=0.9, wspace=0.2, hspace=0.5)

    #plt.tight_layout()
    #plt.show(block=False)

    # titulo da figura
    titulo = "Equilibrium distribution ("+r'$b=$'+str(b)+", "+r'$c=$'+str(c)+\
            ", "+r'$n=$'+str(n)+", "+r'$\beta$='+str(beta)+")"
    plt.suptitle(titulo) 

    # Salva fig
    #nome = "eigenvec_b="+str(b)+"_c="+str(c)+"_alpha"+str(alpha)+"_beta="+str(beta)+".eps"
    nome = "eigenvec_b="+str(b)+"_c="+str(c)+"_alpha"+str(alpha)+"_beta="+str(beta)+".png"
    plt.savefig(nome)

    # Fecha fig
    plt.close()

def funcao(b, delta, n): 
    w0 = 1
    nome_pasta = "b=%.1f" %b
    if not os.path.exists(nome_pasta):
        os.makedirs(nome_pasta)
    for alpha in [0.1, 0.5, 1.0, 2.0]: 
        vbeta = np.arange(0.1, 1.01, 0.1)
        for c in [0.03, 0.15, 0.5, 1., 2., 5.]:
            print "b=%.1f, c=%.2f, alpha=%.1f" %(b,c,alpha)
            nome_arq = "b=%.1f/m_critico_c=%.2f_alpha=%.1f.txt" %(b, c, alpha)
            with open(nome_arq, "w") as f: 
                for beta in vbeta:
                    f.write("{0:.2f}".format(beta))
                    f.write("\t\t")
                    m_c = m_critico(w0,b,c,n,delta,alpha,beta)
                    m_c2 = m_critico2(w0,b,c,n,delta,alpha,beta)
                    f.write("{0:.7f}".format(m_c))
                    f.write("\t\t")
                    f.write("{0:.7f}".format(m_c2))
                    f.write("\n")

## Main
def main():

    w0 = 1.             # fitness basal
    b = 10.              # beneficio gerado por altruista
    c = 2.              # custo para um altruista
    n = 26
    delta = 0.01
    alpha = 0.1

    b = 0.0
    c = 0.03
    beta = 0.1
    #print m_critico(w0,b,c,n,delta,alpha,beta)
    
    #b = 0.
    #for c in [0.03, 0.15, 0.5, 1.0, 2.0, 5.0]:
    #    for beta in np.arange(0.1,1.1,0.1):
    #            autovetores(n,b,c,w0,alpha,beta)
    #autovetores(n,b,c,w0,alpha,beta)
   # vb = [0., 1., 2., 10.]
   # nc = 4
   # Parallel(n_jobs=nc)(delayed(funcao)(b,delta,n) for b in vb)

    #for m in [0.393, 0.176, 0.082]:
    for m in [0.393, 0.176, 0.082, 0.288, 0.103]:
        #encontra_beta(1., 0., 0.03, n, delta, m)
        encontra_alpha(w0, b, c, n, delta, m)
if  __name__ =='__main__':
    main()
