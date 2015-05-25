#coding: utf-8

#############################
# Testando o programa TLFW  #
#############################
import logging
import numpy as np
from retest import gera_simulacao
from simTeorico.main import m_critico, m_critico2

LOG_FILENAME = 'teste.log'
logging.basicConfig(filename=LOG_FILENAME, filemode='a',\
                format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s', \
                datefmt='%H:%M:%S',level=logging.INFO)
logger = logging.getLogger(__name__)

# Como se trata do programa TLFW:
beta = 0.0

# Para todas as sim: N = 5000 (grupos) e n = 26 (indivíduos)
# Além disso, mu = 0.0001, delta = 0.01
N = 5000
n = 26
mu = 0.0001
delta = 0.01
pA = 1.0
pmig = 0.8

# Teste para alpha = 1.0,  b = 0.0 e c variável, conforme cvalores
alpha = 1.
b = 0.
cvalores = [0.03, 0.15, 0.5, 1., 2.]

# Numero de pontos por intervalo [0,1] 
npt = 10

for c in cvalores:
    logger.info(u"....Calculando c = %.2f....." %c)
    m_c = m_critico(1,b,c,n,delta,delta,beta)
    m_c2 = m_critico2(1,b,c,n,delta,delta,beta)
    logger.info(u"A taxa de migração crítica é %.3f (ou %.3f - interpol" %(m_c, m_c2))
    #logger.info(u"Logo: abaixo de mig = %.3f, há chance de emergência de altruísmo. Para taxas de migrações mais altas, a emergência se torna implausível" %m_c)
    vec_m1 = np.arange(0.,m_c,1./npt)
    print vec_m1
    ini2 = 0. if vec_m1.size == 0 else vec_m1[-1]
    vec_m2 = np.arange(ini2+1./npt, 1.+1./npt, 1./npt)[::-1]

    for pmig in vec_m1:
        logger.info(u"Para pmig = %.3f - abaixo do valor de m_c" %pmig)
        n_geracoes = gera_simulacao(N, n, 0., b, c, delta, mu, alpha, beta, pmig)
        if n_geracoes >= 5000:
            logger.info(u"Teste inconclusivo")
        else:
            logger.info(u"O numero de geracoes foi: %d" %n_geracoes)

    for pmig in vec_m2:
        logger.info(u"Para pmig = %.3f - acima do valor de m_c" %pmig)
        n_geracoes = gera_simulacao(N, n, 1., b, c, delta, mu, alpha, beta, pmig)
        if n_geracoes >= 5000:
            logger.info(u"Teste inconclusivo")
        else:
            logger.info(u"O numero de geracoes foi: %d\n" %n_geracoes)

# Teste para alpha = 1.0,  b = 2.0 e c variável, conforme cvalores
alpha = 1.
b = 2.
cvalores = [0.03, 0.15, 0.5, 1., 2.]

# Numero de pontos por intervalo [0,1] 
npt = 10

for c in cvalores:
    logger.info(u"....Calculando c = %.2f....." %c)
    m_c = m_critico(1,b,c,n,delta,delta,beta)
    m_valores.append(m_c)
    logger.info(u"A taxa de migração crítica é %.3f" %m_c)
    #logger.info(u"Logo: abaixo de mig = %.3f, há chance de emergência de altruísmo. Para taxas de migrações mais altas, a emergência se torna implausível" %m_c)
    vec_m1 = np.arange(0.,m_c,1./npt)
    vec_m2 = np.arange(vec_m1[len(vec_m1-1)]+1./npt, 1., 1./npt)[::-1]

    for pmig in vec_m1:
        logger.info(u"Para pmig = %.3f - abaixo do valor de m_c" %pmig)
        n_geracoes = gera_simulacao(N, n, 0., b, c, delta, mu, alpha, beta, pmig)
        if n_geracoes >= 5000:
            logger.info(u"Teste inconclusivo")
        else:
            logger.info(u"O numero de geracoes foi: %d" %n_geracoes)

    for pmig in vec_m2:
        logger.info(u"Para pmig = %.3f - acima do valor de m_c" %pmig)
        n_geracoes = gera_simulacao(N, n, 1., b, c, delta, mu, alpha, beta, pmig)
        if n_geracoes >= 5000:
            logger.info(u"Teste inconclusivo")
        else:
            logger.info(u"O numero de geracoes foi: %d\n" %n_geracoes)

