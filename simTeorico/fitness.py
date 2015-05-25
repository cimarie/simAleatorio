# -*- coding: utf-8 -*-
import numpy as np

# Fitness - recebe o tipo (0 ou 1 - altruista ou nao altruista)
# o numero de altruistas no grupo e o numero total de individuos no grupo
def fitness(tipo,k,n,delta,b,c,w0):
    freqA = np.array(k)/float(n)
    return w0 + delta*(b*freqA-(1-tipo)*c)

# Fitness medio - recebe k e n
def fitness_m(k,n,delta,b,c,w0):
    freq_A = np.array(k/float(n))
    return freq_A*fitness(0,k,n,delta,b,c,w0)+(1-freq_A)*fitness(1,k,n,delta,b,c,w0)
