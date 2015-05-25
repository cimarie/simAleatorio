from selmig import *
import eigenvalue as eg
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import pylab

# Funcao rhoc: calcula (autovalor dominante - 1) a partir dos parametros
def rhoc(c,w0,b,m,n,deltaf,deltac,beta):

    S,M = initSM(n,m,deltaf,deltac,b,c,w0,beta)
    return eg.av_dominante(np.dot(M,S))-1

# Encontra m critico ate o qual todos os autovalores sao maiores que 1 (se 0<m<1, usa o metodo da secante para achar raizes)
def c_critico(w0,b,m,n,deltaf,deltac,beta):

    fa = rhoc(0.,w0,b,m,n,deltaf,deltac,beta)

    fb = rhoc(1.,w0,b,m,n,deltaf,deltac,beta)

    if(fa*fb>0):
        if(fa>0):
            return 1.
        else:
            return 0.

    # brenth encontra o zero da funcao rhoc no intervalo 
    return brentq(rhoc, 0., 10., xtol=0.001, args=(w0,b,m,n,deltaf,deltac,beta))

def encontra_cc(w0,beta,b,m,opcao):

    # Abre arquivo com os tamanhos dos grupos
    with open('ndados.txt', 'r') as f:
        exemplos = [int(x) for x in f.readline().split()]

    contador1 = 0

    vetor_deltaf = np.arange(0.01,1.0,0.01,dtype=float)

    vetor_custoc = np.empty((len(exemplos),len(vetor_deltaf)))

    label = []

    for n in exemplos:

        print(n)

        label.append('n = ' + str(n))

        contador2 = 0

        for deltaf in vetor_deltaf:

            deltac = deltaf

            vetor_custoc[contador1][contador2] = c_critico(w0,b,m,n,deltaf,deltac,beta)

            contador2 = contador2 + 1

        contador1 = contador1 + 1

    w, h = plt.figaspect(1.2)
    plt.figure(figsize=(w,h),dpi=100)

    if opcao==0:

        # Mostra um grafico o custo critico vs. deltaf
        plt.plot(vetor_deltaf, vetor_custoc[0,:],label=label[0])
        for i in range(contador1-1):
            plt.hold(True)
            plt.plot(vetor_deltaf,vetor_custoc[i+1,:],label=label[i+1])

        # Label
        plt.xlabel("deltaf")
        plt.ylabel("c critico")

        # Titulo
        titulo = "Custo critico ("+"b = "+str(b)+", m = "+str(m) + ", beta = " + str(beta) + ")"
        plt.title(titulo)

        plt.autoscale()

        # Limite e intervalos do eixo x
        #pylab.xlim([0.0,1.0])
        #pylab.xticks(np.arange(0.0,1.1,0.1))

        # Limite e intervalos do eixo y
        #pylab.ylim([-0.01,0.5])
        #pylab.yticks(np.arange(0.,0.5,0.05))

        # Legenda
        pylab.legend(loc='upper right')
        plt.hold(False)

        # Grid
        plt.grid(True)

        plt.tight_layout()

        # Salva fig
        nome = "c_versus_deltaf_beta="+str(beta) + "_m="+ str(m)
        plt.savefig(nome+".png")

        # Fecha fig
        plt.close()

    else:
        nome = "c_versus_deltaf_beta="+str(beta) + "_m="+ str(m)
        nome = nome + ".txt"
        arq = open(nome,'w')

        for cont in range(len(vetor_deltaf)):
            arq.write(str(vetor_deltaf[cont]))
            for i in range(len(exemplos)):
                arq.write("\t")
                arq.write(str(vetor_custoc[i,cont]))

            arq.write("\n")

        arq.close()

# Gera figura beta(m) sem as cores ciano magenta
def beta_c(w0, b, m, n, deltaf,deltac,opcao):

    r = 5000

    vetor_beta = np.arange(0.0,1.+1./r,1./r,dtype=float)

    vetor_custo = np.empty(len(vetor_beta))

    for i in range(len(vetor_beta)):

        vetor_custo[i]=c_critico(w0,b,m,n,deltaf,deltac,vetor_beta[i])
        
        #print i

    w, h = plt.figaspect(1)
    plt.figure(figsize=(w,h),dpi=100)

    if opcao==0:

        # Mostra um grafico a migracao critica vs. deltaf
        plt.plot(vetor_beta,vetor_custo)

        # Label
        plt.xlabel("beta")
        plt.ylabel("c_critico")

        # Titulo
        titulo = "Custo critico x beta ("+"b = "+str(b)+", m = "+str(m) + ", n = " + str(n) + ", deltaf = " + str(deltaf) + ")"
        plt.title(titulo)

        plt.autoscale()

        # Limite e intervalos do eixo x
        #pylab.xlim([-0.1,1.1])
        #pylab.xticks(np.arange(0.0,1.05,0.2))

        # Limite e intervalos do eixo y
        #pylab.ylim([-0.1,1.1])
        #pylab.yticks(np.arange(0.0,1.05,0.2))

        # Grid
        plt.grid(True)

        plt.tight_layout()

        # Salva fig
        nome = "c_vs_beta_deltaf="+str(deltaf) + "_m=" + str(m)
        plt.savefig(nome+".png")

        # Fecha fig
        plt.close()

    else:
        nome = "c_vs_beta_deltaf="+str(deltaf) + "_m=" + str(m)
        nome = nome + ".txt"
        arq = open(nome,'w')

        for cont in range(len(vetor_custo)):
            arq.write(str(vetor_beta[cont]))
            arq.write("\t")
            arq.write(str(vetor_custo[cont]))
            arq.write("\n")

        arq.close()


