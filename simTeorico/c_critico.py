from selmig import *
import eigenvalue as eg
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from pylab import *

rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'cm'

# Funcao rhoc: calcula (autovalor dominante - 1) a partir dos parametros
def rhoc(c,w0,b,m,n,delta,alpha,beta):

    S,M = initSM(n,m,delta,alpha,b,c,w0,beta)
    return eg.av_dominante(np.dot(M,S))-1

# Encontra c critico ate o qual todos os autovalores sao maiores que 1 (se 0<m<1, usa o metodo da secante para achar raizes)
def c_critico(w0,b,m,n,delta,alpha,beta):

    fa = rhoc(0.,w0,b,m,n,delta,alpha,beta)

    fb = rhoc(10.,w0,b,m,n,delta,alpha,beta)

    if(fa*fb>0):
        if(fa>0):
            return 10.
        else:
            return 0.

    # brenth encontra o zero da funcao rhoc no intervalo 
    return brentq(rhoc, 0., 10., xtol=0.001, args=(w0,b,m,n,delta,alpha,beta))

#def encontra_cc(w0, alpha, beta, b, m):
#    
#    pass
#

def cc_delta(w0,alpha,beta,b,m,opcao):

    # Abre arquivo com os tamanhos dos grupos
    with open('ndados.txt', 'r') as f:
        exemplos = [int(x) for x in f.readline().split()]

    contador1 = 0

    vetor_delta = np.arange(0.01,1.0,0.01,dtype=float)

    vetor_custoc = np.empty((len(exemplos),len(vetor_delta)))

    label = []

    for n in exemplos:

        print(n)

        label.append(r'$n = $' + str(n))

        contador2 = 0

        for delta in vetor_delta:

            alpha = delta

            vetor_custoc[contador1][contador2] = c_critico(w0,b,m,n,delta,alpha,beta)

            contador2 = contador2 + 1

        contador1 = contador1 + 1

    w, h = plt.figaspect(1.2)
    plt.figure(figsize=(w,h),dpi=100)

    if opcao==0:

        # Mostra um grafico o custo critico vs. delta
        plt.plot(vetor_delta, vetor_custoc[0,:],label=label[0])
        for i in range(contador1-1):
            plt.hold(True)
            plt.plot(vetor_delta,vetor_custoc[i+1,:],label=label[i+1])

        # Label
        plt.xlabel(r"$delta$")
        plt.ylabel(r"$critical c$")

        # Titulo
        titulo = "Critical cost ("+"b = "+str(b)+", m = "+str(m) + ", beta = " + str(beta) + ")"
        plt.title(titulo)

        plt.autoscale()

        # Limite e intervalos do eixo x
        #xlim([0.0,1.0])
        #xticks(np.arange(0.0,1.1,0.1))

        # Limite e intervalos do eixo y
        #ylim([-0.01,0.5])
        #yticks(np.arange(0.,0.5,0.05))

        # Legenda
        legend(loc='upper right')
        plt.hold(False)

        # Grid
        plt.grid(True)

        plt.tight_layout()

        # Salva fig
        nome = "c_versus_delta_beta="+str(beta) + "_m="+ str(m)
        plt.savefig(nome+".png")

        # Fecha fig
        plt.close()

    else:
        nome = "c_versus_delta_beta="+str(beta) + "_m="+ str(m)
        nome = nome + ".txt"
        arq = open(nome,'w')

        for cont in range(len(vetor_delta)):
            arq.write(str(vetor_delta[cont]))
            for i in range(len(exemplos)):
                arq.write("\t")
                arq.write(str(vetor_custoc[i,cont]))

            arq.write("\n")

        arq.close()

# Gera figura beta(m) sem as cores ciano magenta
def beta_c(w0, b, m, n, delta,alpha,opcao=0):

    r = 100

    vetor_beta = np.arange(0.0,1.+1./r,1./r,dtype=float)

    vetor_custo = np.empty(len(vetor_beta))

    for i in range(len(vetor_beta)):

        vetor_custo[i]=c_critico(w0,b,m,n,delta,alpha,vetor_beta[i])
        
        #print i

    w, h = plt.figaspect(1)
    plt.figure(figsize=(w,h),dpi=100)

    if opcao==0:

        # Mostra um grafico a migracao critica vs. delta
        plt.plot(vetor_beta,vetor_custo)

        # Label
        plt.xlabel(r"$\beta$")
        plt.ylabel(r"critical $c$")

        # Titulo
        titulo = "Critical cost vs. " + r"$\beta$" + "("+r"$b = $"+str(b)+r", $m$ = "+str(m) + r", $n$ = " + str(n) + r", $\delta$ = " + str(delta) + ")"
        plt.title(titulo)

        plt.autoscale()

        # Limite e intervalos do eixo x
        #xlim([-0.1,1.1])
        #xticks(np.arange(0.0,1.05,0.2))

        # Limite e intervalos do eixo y
        #ylim([-0.1,1.1])
        #yticks(np.arange(0.0,1.05,0.2))

        # Grid
        plt.grid(True)

        plt.tight_layout()

        # Salva fig
        nome = "c_vs_beta_delta="+str(delta) + "_m=" + str(m)
        plt.savefig(nome+".png")

        # Fecha fig
        plt.close()

    else:
        nome = "c_vs_beta_delta="+str(delta) + "_m=" + str(m)
        nome = nome + ".txt"
        arq = open(nome,'w')

        for cont in range(len(vetor_custo)):
            arq.write(str(vetor_beta[cont]))
            arq.write("\t")
            arq.write(str(vetor_custo[cont]))
            arq.write("\n")

        arq.close()


