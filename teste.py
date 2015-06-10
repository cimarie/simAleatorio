#coding: utf-8

#############################
# Testando o programa TLFW  #
#############################
from joblib import Parallel, delayed
import multiprocessing as mp
import logging
import numpy as np
from retest2 import gera_simulacao
from simTeorico.main import m_critico, m_critico2

LOG_FILENAME = 'teste_cvaloresmaiores.log'
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
fh = logging.FileHandler(LOG_FILENAME)
fh.setLevel(logging.INFO)
#fh.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

################## PARAMETROS GLOBAIS #########################
# Como se trata do programa TLFW:
beta = 0.0
alpha = 1.0  # nao importa

# Para todas as sim: N = 5000 (grupos) e n = 26 (indivíduos)
# Além disso, mu = 0.0001, delta = 0.01
N = 5000
n = 26
mu = 0.0001
delta = 0.01
###############################################################

def bobo(b, c, m_c, pmig):

    pA = 1. if pmig > m_c else 0.

    l_its = [gera_simulacao(N, n, pA, b, c, delta, mu, alpha, beta, pmig) for i in xrange(5)]

    # Mostra quantas simulacoes acabaram antes do maximo e a media de iteracoes que elas levaram
    v = np.array(l_its)<5000
    if np.count_nonzero(v)>0:
        w = np.take(l_its, np.nonzero(v)[0])
        n_geracoes = sum(w)/len(w)
    else:
        n_geracoes = 0 

    return b, c, pmig, np.count_nonzero(v), n_geracoes, l_its[0], l_its[1], l_its[2], l_its[3], l_its[4]

def testa(b, cvalores, numpt):
    for c in cvalores:
        #logger.info(u"....Calculando b = %.2f, c = %.2f....." %(b, c))
        m_c = m_critico(1,b,c,n,delta,delta,beta)
        m_c2 = m_critico2(1,b,c,n,delta,delta,beta)
        logger.info(u"A taxa de migração crítica é %.3f (ou %.3f - interpol)" %(m_c, m_c2))
        #logger.info(u"Logo: abaixo de mig = %.3f, há chance de emergência de altruísmo. Para taxas de migrações mais altas, a emergência se torna implausível" %m_c)

        vec_m = np.arange(0.,1.+1./(10*numpt), 1./numpt)
        for mm in vec_m:
            logger.debug("%.3f esta em vec_m1" %mm)

        logger.debug("\n\n")
        nc = mp.cpu_count()-1
        res1 = Parallel(n_jobs=nc)(delayed(bobo)(b,c,m_c,pmig) for pmig in vec_m)

        l1 = np.array(zip(*res1)).T

        for elem in l1:
            benefit = elem[0]
            cost = elem[1]
            pmig = elem[2]
            palavra = u"abaixo" if pmig<m_c else u"acima"
            logger.info(u"Para b = %.2f, c = %.2f e pmig = %.3f - %s do valor de m_c: %.3f" %(benefit, cost, pmig, palavra, m_c))
            n_fim = elem[3]
            logger.info(u"O numero de simulacoes concluidas é: %d" %n_fim)
            n_geracoes = elem[4]
            if n_geracoes:
                logger.info(u"O numero medio de geracoes foi: %d\n" %n_geracoes)
                logger.info(u"\t%d\t%d\t%d\t%d\t%d\n" %(elem[5],elem[6],elem[7],elem[8],elem[9]))
            else:
                logger.info(u"Teste inconclusivo\n")

def main():
    logger.info(u"Parâmetros fixos nessa simulação: N=%d, n=%d, delta=%.3f,\
            \n\t\tmu=%.4f, alpha=%.1f, beta=%.2f" \
            %(N, n, delta, mu, alpha, beta))

    # Numero de pontos por intervalo [0,1] 
    npt = 10

    # Teste para alpha = 1.0,  b = 0.0 e c variável, conforme cvalores
    benefit = 0.
    #vcost = [0.03, 0.15, 0.5, 1., 2.]
    vcost = [4., 5., 10.]
    testa(benefit,vcost,npt)
    
    # Teste para alpha = 1.0,  b = 2.0 e c variável, conforme cvalores
    benefit = 2.
    testa(benefit,vcost,npt)

    # Teste para alpha = 1.0,  b = 5.0 e c variável, conforme cvalores
    benefit = 5.
    testa(benefit,vcost,npt)

if __name__ == '__main__':
    main()
