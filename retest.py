# coding: utf-8
import logging
import click
import numpy as np
import numpy.random
import random
import time
import scipy.stats
from scipy import sparse
from collections import Counter
from itertools import chain,izip

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#@click.command()
#@click.option('--gnum', default=5000, help='numero de grupos')
#@click.option('--inum', default=26, help='numero de individuos por grupos')
#@click.option('--pa', default=1., help='fracao de altruistas na pop inicial')
#@click.option('--b', default=0., help='beneficio provido pelo altruista em suas interacoes')
#@click.option('--c', default=10., help='custo a que o altruista incorre')
#@click.option('--delta', default=0.01, help='forca de selecao')
#@click.option('--mutacao', default=0.0001, help='taxa de mutacao')
#@click.option('--alpha', default=2., help='prevalencia do altruista em batalhas')
#@click.option('--beta', default=0.0, help='probabilidade de ocorrencia de guerra')
#@click.option('--pmig', default=0.0, help='probabilidade de migracao intergrupos')
def gera_simulacao(gnum, inum, pa, b, c, delta, mutacao, alpha, beta, pmig):

    logger.info(u"Começando a simulação")
    logger.info(u"Parâmetros: N=%d, n=%d, pA=%.2f, b=%.2f, c=%.2f, delta=%.3f,\
        \n\t\tmu=%.4f, alpha=%.1f, beta=%.2f, pmig=%.2f" \
        %(gnum, inum, pa, b, c, delta, mutacao, alpha, beta, pmig))

    A = gnum*inum*pa
    grupos,lfit,lfit_m,mpvencer = initSim(gnum,inum,A,b,c,delta,alpha)
    x = time.time()+random.random()
    return simula(gnum,inum,mutacao,beta,pmig,grupos,lfit,lfit_m,mpvencer,x)

def simula(N, n, PM, beta, pmig, grupos, listafitness, listafitness_m, mpvencer, x):

    random.seed(time.time()+random.random()+x)

    IT = 5002
    precisao = 0.01

    AL = [] 
    AL.append(np.count_nonzero(grupos)/(N*n))
    crit = 0. if AL[0] > (1.-precisao) else 1.

    # Para cada periodo, os grupos entram em conflito e se reproduzem, e
    # os individuos sofrem mutacao e migram entre os grupos
    for it in range(1,IT):
        if abs(AL[it-1]-crit)<precisao:
            print "Acabou na geracao ", it -1
            break
        # 
        knums = [np.count_nonzero(line) for line in grupos]
        glabels = conflito(N,knums,beta,listafitness_m, mpvencer) if N>1 \
                    else knums
        grupos = reproducao_ind(N,n,listafitness,listafitness_m,glabels)
        grupos = mutacao(N,n,PM,grupos)
        grupos = migracao(N,n,grupos,pmig)
        freqA = float(np.count_nonzero(grupos))/(N*n)
        AL.append(freqA)

        #print it, "\t----------->\t", freqA

    return it-1

def initSim(N, n, A, b, c, delta, alpha):    

    grupos = np.zeros((N,n))
    grupos = inicializa(N,n,A,grupos) if A<(N*n) else np.ones((N,n))

    listafitness = cria_listaf(b,c,n,delta)
    listafitness_m = cria_listafm(b,c,n,delta)
    mpvencer = matriz_vencedores(alpha, n)
    
    return grupos, listafitness, listafitness_m, mpvencer

# Inicializa matriz grupos (cada linha é um grupo, com n individuos do tipo 0 ou 1)
def inicializa(N, n, A, grupos):

    while A > 0:
        lin = random.randint(0,N-1)
        col = random.randint(0,n-1)
        if grupos[lin][col]==0:
            grupos[lin][col]=1
            A -= 1

    return sparse.csr_matrix(grupos)
 
def mutacao(N, n, PM, grupos):

    flips = scipy.stats.bernoulli.rvs(PM, size=(N,n))
    return (grupos+flips)%2

def cria_listaf(b, c, n, delta):
    return [fitness(0,n,k,b,c,delta) for k in range(0,int(n)+1)] 

def cria_listafm(b, c, n, delta):
    return [fitness_m(n,k,b,c,delta) for k in range(0,int(n)+1)]

def fitness(tipo, n, k, b, c, delta):
    return 1 + delta*(b*float(k)/n - (1-tipo)*c)

def fitness_m(n, k, b, c, delta):
    freq = float(k)/n
    return freq*fitness(0,n,k,b,c,delta)+(1-freq)*fitness(1,n,k,b,c,delta)

# Reproducao de grupo (só ocorre se não há conflito, beta = 0)
def reproducao_grupo(N, knums, lfitm):

    # Reproducao de grupo
    distk = Counter(knums)

    xk = distk.keys()
    pk = [distk[chave]*lfitm[chave] for chave in xk]
    pk = pk/np.sum(pk)
    distw = scipy.stats.rv_discrete(name='distw', values=(xk,pk))

    glabels = distw.rvs(size=N)

    return glabels

# Reproducao individual (ocorre apos o conflito/reproducao de grupo)
def reproducao_ind(N, n, lfit, lfitm, glabels):

    # Reproducao individual
    aux = np.zeros((N,n))

    for i in range(N):
        k = glabels[i] 
        wA = lfit[k] 
        wm = lfitm[k]
        filhosA = numpy.random.binomial(n, k*wA/(n*wm))
        aux[i][0:filhosA] = 1

    map(numpy.random.shuffle, aux)
    return aux 

# Probabilidade de i vencer j
def pVencer(i, j, alpha, n):
    arg1 = 4.*alpha*i/n
    arg2 = 4.*alpha*j/n
    return np.exp(arg1)/(np.exp(arg1)+np.exp(arg2))

# Constroi matriz com as probabilidades de vencer (linha vencer coluna)
def matriz_vencedores(alpha, n):
    return np.fromfunction(lambda i,j: pVencer(i,j,alpha,n), (n+1,n+1))

# Conflito entre os grupos
def conflito(N, knums, beta, lfitm, mpvencer):

    if beta == 0:
        return reproducao_grupo(N, knums, lfitm)

    else:
        # Calcula numero de grupos que se envolvem em conflitos
        lconts = numpy.random.binomial(N,beta)
        indices = numpy.random.permutation(N)
        lconts = lconts if lconts%2==0 else lconts-1
        lista_conflitos = izip(indices[0:lconts/2], indices[lconts/2:lconts])

        for elem in lista_conflitos:
            i,j = elem
            ki,kj = knums[i], knums[j]
            venc,perd = (i,j) if mpvencer[ki,kj] > random.random() else (j,i)
            knums[perd]=knums[venc]

    return knums   

# Migracao entre os grupos
def migracao(N, n, grupos, pmig):
    
    if pmig > 0.:

        # Lista com numero de migrantes por grupo
        nmigs = numpy.random.binomial(n,pmig,N) 

        # Lista com os migrantes
        migs = list(chain.from_iterable([grupos[i][0:nmigs[i]] for i in range(N)]))
        numpy.random.shuffle(migs)

        cont = 0
        for i in range(N):
            grupos[i][0:nmigs[i]] = migs[cont:cont+nmigs[i]]
            cont += nmigs[i]

        map(numpy.random.shuffle, grupos)

    return grupos

#if __name__ == "__main__":
#    gera_simulacao()
