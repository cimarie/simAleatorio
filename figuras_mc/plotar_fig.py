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
def imprime_figura(b, c, delta, alpha, n):

    titulo = r"Critical migration rate comparison (model vs. simulation) "
    titulo = titulo + "\n" + r"$b$ = %.2f, $c$ = %.2f, $\delta$ = %.3f, $\alpha$ = %.1f and $n$ = %d" \
            %(b, c, delta, alpha, n)        

    lmarkers = ["+", "*"]
    data2 = []
    labels = []
    nome = "teste_mc_c=%.2f" % c
    m = np.loadtxt(nome+'.txt') 
    df = pd.DataFrame(m, columns=["beta", "mc_sim", "mc_teo"])
    data1 = df["beta"].values
    col_mc1 = df["mc_sim"].values
    col_mc2 = df["mc_teo"].values

    fig = plt.figure(figsize=(7,6), dpi=300)
    ax = fig.add_subplot(111)
    data1 = [float(i) for i in data1]
    col_mc1 = [float(i) for i in col_mc1]
    col_mc2 = [float(i) for i in col_mc2]

    ax.plot(data1, col_mc1, marker=lmarkers[0], label=r"simulation")
    ax.plot(data1, col_mc2, marker=lmarkers[1], label=r"theoretical model")
    plt.title(titulo)
        
   # w, h = plt.figaspect(1)
    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$m_c$")
    plt.xlim([-0.01,1.01])
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
    nome_fig = "comparacao_c=%.2f_delta=%.2f_b=%.1f_alpha=%.1f.png" \
                    %(c,delta,b,alpha)
    plt.savefig(nome_fig, additional_artists = lgd, bbox_inches="tight")

    plt.close()

def main():

    delta = 0.01
    n = 26
    b = 1.0
    alpha = 0.1
    
    for c in [0.03, 0.15, 0.5, 1.0, 2.0, 5.]:
        imprime_figura(b, c, delta, alpha, n)

if __name__ == "__main__": main()
#    cli()
