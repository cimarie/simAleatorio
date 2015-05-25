import numpy as np
import matplotlib.pyplot as plt
import pylab

def devolve_lambda(mu, Fst, b, c, n):

    return (c-b/n)/(4*mu*((Fst/(1-Fst))+(1./n)))

def main():

    c = 3
    b = 5 
    Fst_vec = [0.022,0.075,0.170]
    n = 26 
    label = []
    mu_vec = np.arange(0.0001,0.26,0.00001, dtype=float)
    val_lambda = np.empty((len(Fst_vec), len(mu_vec)))

    for i in range(len(Fst_vec)):
   
        label.append("Fst = " + str(Fst_vec[i]))
    
        for j in range(len(mu_vec)):
            val_lambda[i][j]=devolve_lambda(mu_vec[j],Fst_vec[i],b,c,n)

    plt.plot(mu_vec,val_lambda[0,:],label=label[0])
    plt.hold(True)
    plt.plot(mu_vec,val_lambda[1,:],label=label[1])
    plt.hold(True)
    plt.plot(mu_vec,val_lambda[2,:],label=label[2])

    plt.xlabel(r'$\delta$')
    plt.ylabel(r'$\lambda_A$')

    plt.grid(True)

    pylab.legend()
    plt.hold(False)

    pylab.xlim([0.0,0.05])
    #pylab.xticks(np.arange(0.0,0.26,0.05))

    #pylab.ylim([0.0,3.0])
    #pylab.yticks(np.arange(0.0,3.1,0.5))

#    plt.autoscale(True)

    plt.tight_layout()

    nome = "figura_bowles_c.png"
    plt.savefig(nome)

    plt.close()

if __name__ == '__main__':
    main()
