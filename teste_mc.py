#coding: utf-8

################################
# Testando a simulacao Bowles  #
################################
from joblib import Parallel, delayed
import multiprocessing as mp
import logging
import numpy as np
from retest2 import gera_simulacao
from simTeorico.main import m_critico, m_critico2

################## PARAMETROS GLOBAIS #########################
alpha = 0.1 

# Para todas as sim: N = 5000 (grupos) e n = 26 (indivíduos)
# Além disso, mu = 0.0001, delta = 0.01
N = 5000
n = 26
mu = 0.0001
delta = 0.01
max_it = 5
###############################################################

LOG_FILENAME = 'teste_mc.log'
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
fh = logging.FileHandler(LOG_FILENAME)
#fh.setLevel(logging.INFO)
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)
logger.setLevel(logging.DEBUG)

def sim(b, c, beta, pA, pmig, i):
    return gera_simulacao(N, n, pA, b, c, delta, mu, alpha, beta, pmig, max_it)

def bobo(b, c, beta, m_c, pmig):

    pA = 0.
    #pA = 1. if pmig > m_c else 0.
    nc = mp.cpu_count() -1 
    l_its = Parallel(n_jobs=nc)(delayed(sim)(b,c,beta,pA,pmig,i) for i in xrange(10))

    # Mostra quantas simulacoes acabaram antes do maximo e a media de iteracoes que elas levaram
    v = np.array(l_its)<max_it
    if np.count_nonzero(v)>0:
        w = np.take(l_its, np.nonzero(v)[0])
        n_geracoes = sum(w)/len(w)
    else:
        n_geracoes = max_it 
    return np.count_nonzero(v), n_geracoes


def testa(b, c, vbeta, numpt):
    nome_arq = "teste_mc_c=%.2f.txt" %c 
    f = open(nome_arq, "w")
    precisao = 0.0001

    for beta in vbeta:
        m_c_pred = m_critico(1,b,c,n,delta,alpha,beta)
        #logger.info(u"Logo: abaixo de mig = %.3f, há chance de emergência de altruísmo. Para taxas de migrações mais altas, a emergência se torna implausível" %m_c)
        print "BETA: %.f" %beta
        continua = False
        primeiro, ultimo = 0., 1.
        while (continua) or (ultimo-primeiro>precisao):
            logger.debug("primeiro: %f" %primeiro)
            logger.debug("ultimo: %f" %ultimo)
            vec_m = np.arange(primeiro,ultimo+1./(10*numpt), 1./numpt)
            for i in xrange(len(vec_m)):
                pmig = vec_m[i]
                print "i: %d, pmig: %f" %(i,pmig)
                logger.debug("Migracao: %f" %pmig)
                nfim, ngeracoes = bobo(b,c,beta,m_c_pred,pmig)
                logger.debug("Numero de simulacoes ate o fim: %d" %nfim)
                if nfim == 0:
                    if i == 0:
                        continua=False
                    else:
                        primeiro, ultimo = vec_m[i-1], vec_m[i]
                        numpt = 10*numpt
                    break

        m_c_real = pmig
        f.write(str(beta))
        f.write("\t")
        f.write(str(m_c_real))            
        f.write("\t")
        f.write(str(m_c_pred))
        f.write("\n")

    f.close()

def main():

    benefit = 0.
    logger.info(u"Parâmetros fixos nessa simulação: N=%d, n=%d, b=%.2f, \
         delta=%.3f, \n\t\tmu=%.4f, alpha=%.1f" \
            %(N, n, benefit, delta, mu, alpha))

    # Numero de pontos por intervalo [0,1] 
    npt = 10

    vbeta = np.arange(0.,1.01, 0.1)

    # Teste para b = 0.0 e c = 0.5 
    vcost = [0.03, 0.15, 0.5, 1.0, 2.0, 5.0]
    for cost in vcost:
        testa(benefit,cost,vbeta,npt)

if __name__ == '__main__':
    main()
