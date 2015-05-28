#coding: utf-8

#############################
# Testando o programa TLFW  #
#############################
from joblib import Parallel, delayed
import logging
import numpy as np
from retest2 import gera_simulacao
from simTeorico.main import m_critico, m_critico2

LOG_FILENAME = 'teste.log'
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
fh = logging.FileHandler(LOG_FILENAME)
fh.setLevel(logging.INFO)
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
pA = 1.0
pmig = 0.8
###############################################################

def bobo(pA, b, c, pmig):

    l_its = [gera_simulacao(N, n, pA, b, c, delta, mu, alpha, beta, pmig) for i in xrange(5)]

    # Mostra quantas simulacoes acabaram antes do maximo e a media de iteracoes que elas levaram
    v = np.array(l_its)<5001
    if np.count_nonzero(v)>0:
        w = np.take(l_its, np.nonzero(v)[0])
        n_geracoes = sum(w)/len(w)
    else:
        n_geracoes = 0 

    return b, c, pmig, np.count_nonzero(v), n_geracoes

def testa(b, cvalores, numpt):
    for c in cvalores:
        #logger.info(u"....Calculando b = %.2f, c = %.2f....." %(b, c))
        m_c = m_critico(1,b,c,n,delta,delta,beta)
        m_c2 = m_critico2(1,b,c,n,delta,delta,beta)
        #logger.info(u"A taxa de migração crítica é %.3f (ou %.3f - interpol" %(m_c, m_c2))
        #logger.info(u"Logo: abaixo de mig = %.3f, há chance de emergência de altruísmo. Para taxas de migrações mais altas, a emergência se torna implausível" %m_c)
        vec_m1 = np.arange(0.,m_c,1./npt)
        ini2 = 0. if vec_m1.size == 0 else vec_m1[-1]
        vec_m2 = np.arange(ini2+1./npt, 1.+1./(10*npt), 1./npt)[::-1]

        nc = 6
        res1 = Parallel(n_jobs=nc)(delayed(bobo)(0.,b,c,pmig) for pmig in vec_m1)
        res2 = Parallel(n_jobs=nc)(delayed(bobo)(1.,b,c,pmig) for pmig in vec_m2)

        l1 = np.array(zip(*res1)).T
        l2 = np.array(zip(*res2)).T

        for elem in l1:
            benefit = elem[0]
            cost = elem[1]
            pmig = elem[2]
            logger.info(u"Para b = %.2f, c = %.2f e pmig = %.3f - abaixo do valor de m_c:" %(benefit, cost, pmig))
            n_fim = elem[3]
            logger.info(u"O numero de simulacoes concluidas é: %d" %n_fim)
            n_geracoes = elem[4]
            if n_geracoes:
                logger.info(u"O numero medio de geracoes foi: %d\n" %n_geracoes)
            else:
                logger.info(u"Teste inconclusivo\n")

        for elem in l2:
            benefit = elem[0]
            cost = elem[1]
            pmig = elem[2]
            logger.info(u"Para b = %.2f, c = %.2f e pmig = %.3f - acima do valor de m_c:" %(benefit, cost, pmig))
            n_fim = elem[3]
            logger.info(u"O numero de simulacoes concluidas é: %d" %n_fim)
            n_geracoes = elem[4]
            if n_geracoes:
                logger.info(u"O numero medio de geracoes foi: %d\n" %n_geracoes)
            else:
                logger.info(u"Teste inconclusivo\n")

def main():

    logger.info(u"Parâmetros fixos nessa simulação: N=%d, n=%d, pA=%.2f, delta=%.3f,\
            \n\t\tmu=%.4f, alpha=%.1f, beta=%.2f" \
            %(N, n, pA, delta, mu, alpha, beta))
    
    # Teste para alpha = 1.0,  b = 0.0 e c variável, conforme cvalores
    benefit = 0.
    vcost = [0.03, 0.15, 0.5, 1., 2.]
    # Numero de pontos por intervalo [0,1] 
    npt = 10

    testa(benefit,vcost,npt)

    # Teste para alpha = 1.0,  b = 2.0 e c variável, conforme cvalores
    b = 2.
    cvalores = [0.03, 0.15, 0.5, 1., 2.]

    testa(benefit,vcost,npt)

if __name__ == '__main__':
    main()
