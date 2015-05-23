# -*- coding: utf-8 -*-
#eigenvalue.py

import scipy.linalg
import numpy as np

def av_dominante(Matriz):

    Matriz = np.delete(Matriz,0,0)

    Matriz = np.delete(Matriz,0,1)

    av = scipy.linalg.eigvals(Matriz)

    x = np.argsort(-abs(av))

    return np.real(av[x[0]])
