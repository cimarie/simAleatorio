#coding: utf-8

#############################
# Testando o programa TLFW  #
#############################
from joblib import Parallel, delayed
import logging
import numpy as np
from retest2 import gera_simulacao
from simTeorico.main import m_critico, m_critico2

LOG_FILENAME = 'teste_paralelo.log'
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
fh = logging.FileHandler(LOG_FILENAME)
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

############## PARAMETROS GLOBAIS ##################
# Como se trata do programa TLFW:
alpha = 1.0 # Não importa
beta = 0.0

# Para todas as sim: N = 5000 (grupos) e n = 26 (indivíduos)
# Além disso, mu = 0.0001, delta = 0.01
N = 5000
n = 26
mu = 0.0001
delta = 0.01
pA = 1.0
####################################################

def sim(b, c, pmig, i):
    return gera_simulacao(N, n, pA, b, c, delta, mu, alpha, beta, pmig)

def bobo(pA, b, c, pmig):
    l_its = Parallel(n_jobs=5)(delayed(sim)(b,c,pmig,i) for i in xrange(10))

    # Mostra quantas simulacoes acabaram antes do maximo e a media de iteracoes que elas levaram
    v = np.array(l_its)<5001
    if np.count_nonzero(v)>0:
        w = np.take(l_its, np.nonzero(v)[0])
        n_geracoes = sum(w)/len(w)
    else:
        n_geracoes = 0 
    return n_geracoes

def testa(b, cvalores, npt):
    for c in cvalores:
        logger.info(u"....Calculando b = %.2f, c = %.2f....." %(b, c))
        m_c = m_critico(1,b,c,n,delta,delta,beta)
        m_c2 = m_critico2(1,b,c,n,delta,delta,beta)
        logger.info(u"A taxa de migração crítica é %.3f (ou %.3f - interpol)" %(m_c, m_c2))
        #logger.info(u"Logo: abaixo de mig = %.3f, há chance de emergência de altruísmo. Para taxas de migrações mais altas, a emergência se torna implausível" %m_c)
        vec_m1 = np.arange(0.,m_c,1./npt)
        ini2 = 0. if vec_m1.size == 0 else vec_m1[-1]
        vec_m2 = np.arange(ini2+1./npt, 1.+1./(10*npt), 1./npt)[::-1]

        for pmig in vec_m1:
            n_geracoes = bobo(0.,b,c,pmig)
            if n_geracoes == 0:
                logger.info(u"Teste inconclusivo")
                break
            else:
                logger.info(u"O numero de geracoes foi: %d" %n_geracoes)
        for pmig in vec_m2:
            n_geracoes = bobo(1.,b,c,pmig)
            if n_geracoes == 0:
                logger.info(u"Teste inconclusivo")
                break
            else:
                logger.info(u"O numero de geracoes foi: %d" %n_geracoes)

def main():

    # Teste para alpha = 1.0,  b = 2.0 e c variável, conforme cvalores
    benefit = 2.
    vcost = [0.03, 0.15, 0.5, 1., 2.][::-1]

    # Numero de pontos por intervalo [0,1] 
    numpt = 10

    logger.info(u"Parâmetros fixos nessa simulação: N=%d, n=%d, pA=%.2f, delta=%.3f,\
            \n\t\t mu=%.4f, alpha=%.1f, beta=%.2f"\
            %(N, n, pA, delta, mu, alpha, beta))

    testa(benefit, vcost, numpt)

# Teste para alpha = 1.0,  b = 0.0 e c variável, conforme cvalores
    benefit = 0.
    vcost = [0.03, 0.15, 0.5, 1., 2.][::-1]

    testa(benefit, vcost, numpt)

if __name__ == '__main__':
    main()
