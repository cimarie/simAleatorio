#coding: utf-8
# $ python workflow_teste.py RodaTodosTestes --benefit 0.0 --alpha 0.1 --nindiv 26 --ngroup 5000 --mu 0.0001 --delta 0.001

################################
# Testando a simulacao Bowles  #
################################
from joblib import Parallel, delayed
import multiprocessing as mp
import logging
import luigi
import numpy as np
from retest2 import gera_simulacao
from simTeorico.main import m_critico, m_critico2

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

max_it = 5000

def sim(N, n, b, c, delta, mu, alpha, beta, pA, pmig, i):
    return gera_simulacao(N, n, pA, b, c, delta, mu, alpha, beta, pmig, max_it)

def bobo(N, n, b, c, delta, mu, alpha, beta, pmig):

    pA = 0.0
    it = gera_simulacao(N,n,pA,b,c,delta,mu,alpha,beta,pmig,max_it)
    if it < max_it:
        return 1, it 
    it = gera_simulacao(N,n,1.0,b,c,delta,mu,alpha,beta,pmig,max_it)
    if it < max_it:
        logger.info("pA=1 => 0")
        return 0, max_it
    nc = mp.cpu_count() -1 
    l_its = Parallel(n_jobs=nc)(delayed(sim)(N, n, b, c, delta, mu, alpha, beta, pA, pmig, i)\
                                 for i in xrange(10))

    # Mostra quantas simulacoes acabaram antes do maximo e a media de iteracoes que elas levaram
    v = np.array(l_its)<max_it
    if np.count_nonzero(v)>0:
        w = np.take(l_its, np.nonzero(v)[0])
        n_geracoes = sum(w)/len(w)
    else:
        n_geracoes = max_it 
    return np.count_nonzero(v), n_geracoes


class RodaTodosTestes(luigi.Task):

    benefit = luigi.Parameter(default=0.0)
    alpha = luigi.Parameter(default=0.1)
    ngroup = luigi.Parameter(default=5000)
    nindiv = luigi.Parameter(default=26)
    mu = luigi.Parameter(default=0.0001)
    delta = luigi.Parameter(default=0.01)
    
    def requires(self):
        return None
    
    def run(self):
        logger.info(u"Parâmetros fixos nessa simulação: N=%s, n=%s b=%s, \
     delta=%s, \n\t\tmu=%s, alpha=%s" \
        %(self.ngroup, self.nindiv, self.benefit, self.delta, self.mu, self.alpha))
        for c in ["0.03", "0.15", "0.5", "1.0", "2.0", "5.0"]:
            yield RodaTeste(self.benefit, self.alpha, self.ngroup, self.nindiv, self.mu, self.delta, c)
   # def complete(self):
   #     return False

class RodaTeste(luigi.Task):
    
    benefit = luigi.Parameter(default=0.0)
    alpha = luigi.Parameter(default=0.1)
    ngroup = luigi.Parameter(default=5000)
    nindiv = luigi.Parameter(default=26)
    mu = luigi.Parameter(default=0.0001)
    delta = luigi.Parameter(default=0.01)
    cost = luigi.Parameter(default=0.03)

    def requires(self):
        return None

    def run(self):
        alpha = float(self.alpha)
        n = int(self.nindiv)
        N = int(self.ngroup)
        mu = float(self.mu)
        delta = float(self.delta)
        b = float(self.benefit)
        c = float(self.cost)
        precisao = 0.001
        vbeta = np.arange(0.,1.01, 0.1)

        f = self.output().open('w')         

        for beta in vbeta:
            numpt = 10
            m_c_pred = m_critico(1,b,c,n,delta,alpha,beta)
            logger.info("C: %f \t BETA: %f" %(c,beta))
            continua = True 
            reverso = False
            primeiro, ultimo = 0., 1.
            pmig = primeiro
            while (continua) and (ultimo-primeiro>precisao):
                logger.debug("primeiro: %f" %primeiro)
                logger.debug("ultimo: %f" %ultimo)
                vec_m = np.arange(primeiro,ultimo+1./(10*numpt), 1./numpt)
                if reverso:
                    vec_m = vec_m[::-1]
                for i in xrange(len(vec_m)):
                    pmig = vec_m[i]
                    logger.debug("i: %d, pmig: %f" %(i,pmig))
                    nfim, ngeracoes = bobo(N, n, b, c, delta, mu, alpha, beta, pmig)
                    logger.debug("Numero de simulacoes ate o fim: %d" %nfim)
                    if nfim == 0:
                        if i == 0:
                            continua=False
                            logger.debug("Continua modificado: False")
                        else:
                            numpt = 10*numpt
                            primeiro = vec_m[i-1]+1./numpt
                            ultimo = vec_m[i]-1./numpt
                            reverso=False
                        break
                    elif vec_m[i] == max(vec_m):
                        pmig = vec_m[i] 
                        continua=False
                        reverso=True

            m_c_real = pmig
            logger.info("m_c encontrado: %f" %m_c_real)
            f.write(str(beta))
            f.write("\t")
            f.write(str(m_c_real))            
            f.write("\t")
            f.write(str(m_c_pred))
            f.write("\n")
        f.close()

    def output(self):
        return luigi.LocalTarget('teste_mc_c=%.2f.txt' %float(self.cost))

if __name__ == '__main__':
    luigi.run()
