# -*- coding: utf-8 -*-
import click
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import *

rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'cm'

#@click.group()
#def cli():
#    pass
#
#@cli.command()
#@click.option('--b', default=0., help='beneficio provido pelo altruista em suas interacoes; default: 0.')
#@click.option('--delta', default=0.01, help='forca de selecao; default: 0.01')
#@click.option('--alpha', default=1., help='prevalencia do altruista em batalhas; default: 2.')
#@click.option('--n', default=26, help='numero de individuos por grupos; default: 26')
def imprime_figura(b, delta, alpha, n):

    titulo = r"Critical migration rate  "
    titulo = titulo + r"$b$ = %.3f, $\delta$ = %.3f, $\alpha$ = %.3f and $n$ = %d" \
            %(b, delta, alpha, n)        

    lmarkers = ["v", "^", "o", "h", "+", "*"]
    data2 = []
    labels = []
    for c in [0.03, 0.15, 0.5, 1., 2., 5.]:
        nome = "delta=%.2f/b=%.1f/m_critico_c=%.2f_alpha=%.1f" \
                    % (delta, b, c, alpha)
        labels.append(nome[-16:-10])
        m = np.loadtxt(nome+'.txt') 
        df = pd.DataFrame(m, columns=["beta", "mc1", "mc2"])
        data1 = df["beta"].values
        col_mc = df["mc1"].values
        data2 = np.append(data2, col_mc)
    print labels

    fig = plt.figure(figsize=(7,7), dpi=150)
    ax = fig.add_subplot(111)
    data1 = [float(i) for i in data1]
    data2 = np.array([float(i) for i in data2])
    data2 = data2.reshape(6, len(data1))

    for i in xrange(len(labels)):
        ax.plot(data1,data2[i,:],marker=lmarkers[i], label=labels[i])
    plt.title(titulo)
        
   # w, h = plt.figaspect(1)
    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$m_c$")
    plt.xlim([0.,1.01])
   # plt.xticks(np.arange(0.0,1.05,0.2))
    plt.ylim([-0.01,1.01])
   # plt.yticks(np.arange(0.0,1.05,0.2))
   # plt.show()

    # Shrink current axis's height by 10% on the bottom
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #             box.width, box.height * 0.9])
    #pylab.legend(loc='upper right')
    lgd = ax.legend(loc="center", bbox_to_anchor=(0.5,-0.13), fancybox=True, ncol=3)
    #plt.legend()
    plt.grid(True)
    plt.tight_layout()
    nome_fig = "delta=%.2f/m_critico_b=%.1f_alpha=%.1f.png" %(delta,b,alpha)
    plt.savefig(nome_fig, additional_artists = lgd, bbox_inches="tight")

    plt.close()

def main():

    delta = 0.01
    n = 26
    for b in [0., 1.0, 2.0, 10.]:
        for alpha in [0.1, 0.5, 1., 2.]:
            imprime_figura(b, delta, alpha, n)

if __name__ == "__main__": main()
#    cli()
