using Distributions
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport scipy.optimize as so
@pyimport scipy.misc as sm
@pyimport scipy.linalg as sl
@pyimport smtplib

# Main
function main()

    n = 20 #NUMBER OF MEMBERS IN EACH DEME
    N = 5000 #NUMBER OF DEMES
    beta = 0.5 #THE PROBABILITY OF A DEME TO ENTER IN A CONTEST
    LA = 0.2 #PARAMETER THAT INFLUENCES THE PROBABILITY OF WINNING (L)
    b = .0 #BENEFIT FROM AN INTERACTION
    c = 0.6  #COST TO THE ALTRUIST IN THE INTERACTION
    delta = 0.01 #SELECTION STRENGTH
    PM = 0.01 #MUTATION PROBABILITY

    it = Int32[] 

    npt = 10

    deltac=delta
    deltaf=delta

    m_c = m_critico(1,b,c,n,deltaf,deltac,beta)
	
    vec_m1 = [0.:1./npt:m_c]

    vec_m2 = [vec_m1[length(vec_m1)]+1./npt:1./npt:1.]

    vec_m2 = reverse(vec_m2)

    arq6 = open("terminou", "w")
    arq7 = open("naoterminou", "w")

    for PMIG in vec_m1
        A = 0
        # Starting from A=0, we will see if this fixes at 1
        res = Simulate10(N,n,A,LA,b,c,beta,delta,PM,PMIG,0)
        arq =(res[1])? arq6 : arq7
        write(arq, string(PMIG))
        write(arq, "\t")
        write(arq, string(A/(N*n)))
        for elem in res[2]
            write(arq, "\t")
            write(arq, string(elem))
        end
        write(arq, "\n")
    end

    for PMIG in vec_m2
        A = N*n 
        # Starting from A=N*n, we will see if this fixes at 0
        res = Simulate10(N,n,A,LA,b,c,beta,delta,PM,PMIG,1)
        arq =(res[1])? arq6 : arq7
        write(arq, string(PMIG))
        write(arq, "\t")
        write(arq, string(A/(N*n)))
        for elem in res[2]
            write(arq, "\t")
            write(arq, string(elem))
        end
        write(arq, "\n")          
    end

    close(arq6)
    close(arq7)

    texto =  "simulacao terminou\n\n)"

    manda_email(texto)

end

# Run ten simulations and returns the result of the majority of runs
function Simulate10(N,n,A,LA,b,c,beta,delta,PM,PMIG,op)

    contV = 0

    num_its = Int32[]

    for i in 1:11
        res = Simulate(N,n,A,LA,b,c,beta,delta,PM,PMIG,op)
        push!(num_its, res[2])
        
        if(res[1])
            contV+=1
        end
    end

    #= If tied, runs one more simulation to decide
    if contV == 5
        outro = Simulate(N,n,A,LA,b,c,beta,delta,PM,PMIG,op)
        return outro[1], num_its
    # return true if contV is larger than 5 and false, otherwise
    else
        return (contV>5), num_its
    end
    =#

    return (contV>5), num_its

end

# Manda email
function manda_email(mensagem)

    smtpserver = "smtp.live.com"
    AUTHREQUIRED = 1
    smtpuser = "cinthia_marie@hotmail.com"
    smtppass = "kurniawan"

    RECIPIENTS = "cimarie@gmail.com"
    SENDER = "cinthia_marie@hotmail.com"
    s = mensagem

    server = pycall(smtplib.SMTP, PyAny, smtpserver,587)
    pycall(server["ehlo"], PyAny)
    pycall(server["starttls"], PyAny)
    pycall(server["ehlo"], PyAny)
    pycall(server["login"], PyAny, smtpuser,smtppass)
    pycall(server["set_debuglevel"], PyAny, 1)
    pycall(server["sendmail"], PyAny, SENDER,[RECIPIENTS],s)
    pycall(server["quit"], PyAny)

end


# Run one simulation
function Simulate(N,n,A,LA,b,c,beta,delta,PM,PMIG,op)

    IT = 5000
    precisao = 0.001

    #set counter zero
    it = 0

    #GENERATES THE DEMES MATRIX WITH THE GIVEN PARAMETERS OF N and n:
    DEMES = zeros(Int32,N,n)

    #println(DEMES)

    DEMES = Initialize(N,n,A,DEMES)

    #println(DEMES)

    #ALTRUISTS COUNTER:
    AL = Float64[]
    push!(AL,(sum(DEMES))/(N*n))

    crit = false 

    while it < IT && !crit
    
        DEMES = Make_contest(DEMES, N, n, LA, b, c, beta, delta)
        DEMES = Insert_mutation(N, n, PM, DEMES)
        DEMES = Make_migration(N, n, PMIG, DEMES)

        #ALTRUISTS COUNTER:
        push!(AL,(sum(DEMES))/(N*n))

        it += 1

        crit = abs(AL[it+1]-(1-op))<precisao
    end
   
    if crit
        println(it+1)
    end
    
    #ARRAY WITH THE PERIODS OF TIME FOR WHICH # ALTRUISMS WAS MEASURED
    IC = [0:it]

    #println( it+1)

    #println(sum(DEMES)/(N*n))

    #SAVE RESULTS TO A FILE
    Save_file(IC, AL,"sim.txt")

    # Save results to a figure
    #Save_fig(N, n, A, IC, AL, b, beta, PMIG)
    
    return crit, it 
end

function Regression(N,n,A,LA,b,c,beta,delta,PM)

    it = Int32[] 

	vec_m = [0:0.1:1]

	for PMIG in vec_m
	    println(PMIG)
    	push!(it, Simulate(N,n,A,LA,b,c,beta,delta,PM,PMIG))
	end	

    inv_it = [1/it[i] for i in 1:length(it)]

    # do a linear regression
    a,b = linreg(vec_m, inv_it)

    # do a scatter plot
    scatter(vec_m,inv_it)
    plot(vec_m,[a+b*i for i in vec_m])

    savefig("figura.png")

    close()

    arq = open("1.txt","w")

    # keep the values of number of iterations for each m
    for cont in 1:length(vec_m)
        write(arq,string(vec_m[cont]))
        write(arq,"\t")
        write(arq,string(it[cont]))
        write(arq,"\n")
    end

    close(arq)
end

# Initialize the demes
function Initialize(n_groups, n_ind, n_A, DEMES)

    #println("aloo")

    #PUTS THE ALTRUISTS IN THE DEMES MATRIX:
    Altr=n_A
    while Altr>0
        ll=rand(1:n_groups)
        cc=rand(1:n_ind)

        if DEMES[ll,cc] == 0
            DEMES[ll,cc] = 1
            Altr -= 1
        end
    end

    return DEMES
end

# Make the contest between two demes
function Make_contest(DEMES, N, n, LA, b, c, beta, delta)

    #INTERACTIONS BETWEEN EACH DEME:
    D_PERM = shuffle([1:N]) #ARRAY DE DEMES PERMUTADOS
    
    for cont_D in 1:N
        
        D = DEMES[cont_D,:] #DEME D (focus)
        cont_d = D_PERM[cont_D]
        d = DEMES[cont_d,:] #DEME d (hostile)

        #ANALYSES THE POSSIBILITY OF THE CONTEST:
        ##CALCULATE L:
        xD = sum(D)/n # Frequency OF ALTRUISTS IN GROUP D
        xd = sum(d)/n # Frequency OF ALTRUISTS IN GROUP d
        #P=LA*(xD-xd)
        L = (exp(4*xD))/(exp(4*xD)+exp(4*xd))
        ##CALCULATE PA, THE PROBABILITY OF AN INDIVIDUAL (D or d) BE AN ALTRUIST:
        W_AD = 1 + delta*(-c + b*xD) #FITNESS OF D ALTRUISTS
        W_Ad = 1 + delta*(-c + b*xd) #FITNESS OF d ALTRUISTS
        W_ND = 1 + delta*b*(1-xD) #FITNESS OF D NON-ALTRUISTS
        W_Nd = 1 + delta*b*(1-xd) #FITNESS OF d NON-ALTRUISTS
        W_mD = W_AD*xD + W_ND*(1-xD) #AVERAGE FITNESS OF D
        W_md = W_Ad*xd + W_Nd*(1-xd) #AVERAGE FITNESS OF d

        PAD = (xD*W_AD)/W_mD
        PAd = (xd*W_Ad)/W_md

        ## THE CONTEST HAPPENS: WINNER TAKES 2 DEMES, LOSER TAKES 0 DEMES
        if rand() < beta
            
            if rand() < L
                # Deme D wins
                winRate = PAD
            else
                # Deme d wins
                winRate = PAd
            end
        
            ## CALCULATE THE OFFSPRING of the winner:
            OffD = rand(Binomial(n, winRate))
            Offd = rand(Binomial(n, winRate))
            
            ## THE CONTEST DOESN'T HAPPEN: each group reproduces in its own deme
        else
            ## CALCULATE THE OFFSPRING OF D:
            OffD = rand(Binomial(n, PAD))
            ## CALCULATE THE OFFSPRING OF d:
            Offd = rand(Binomial(n, PAd))
        
        end

        ## Reproduction in D
        DEMES[cont_D,:] = vcat(fill(1,OffD),fill(0,n-OffD))
        DEMES[cont_D,:] = shuffle(vec(DEMES[cont_D,:]))
        ## Reproduction in d
        DEMES[cont_d,:] = vcat(fill(1,Offd),fill(0,n-Offd))
        DEMES[cont_d,:] = shuffle(vec(DEMES[cont_d,:]))

    end

    return DEMES
end

function Insert_mutation(n_groups, n_ind, mu, DEMES)

    #MUTATION
    flips = rand!(Bernoulli(mu),Array(Int32,n_groups,n_ind))
    DEMES = mod(DEMES+flips,2)

    return DEMES
end

function Make_migration(N, n, PMIG, DEMES)
    #MIGRATION
    VDEME = Int32[]
    VPOS = Int32[]
    VMIG = Int32[]
    ##CHOOSES RANDOMLY SOME INDIVIDUALS FROM EACH DEME TO MIGRATE AND PUTS IN VMIG:
    for group in 1:N
        
        migrators = rand(Binomial(n, PMIG))
        migrated = 0
        
        while migrators > 0
            
            migrated += 1

            # Take the first ones to migrate
            push!(VMIG, DEMES[group,migrated])
            # Fill VDEME with the groups that have blank spaces
            push!(VDEME,group)
            migrators -= 1
        end
        
        if migrated > 0
            VPOS=vcat(VPOS,[1:migrated])
        end

    end        
    
    if length(VMIG) > 0
      
        ##SHUFFLE VMIG
        shuffle(VMIG)
        
        ##REDISTRIBUTE THE INDIVIDUALS INTO THE DEMES, IN ORDER TO MANTAIN ITS SIZE n
        
        for cont in 1:length(VMIG)
            DEMES[VDEME[cont],VPOS[cont]]=VMIG[cont]
        end

        for line in 1:N
            DEMES[line,:] = shuffle(vec(DEMES[line,:]))
        end
    
    end

    return DEMES
end

function Save_fig(N, n, A, IC, AL, b, beta, PMIG)

    plot(IC, AL)
    title("FREQ. OF A ( #A = $A, b=$b, beta=$beta, m = $PMIG )")
    xlabel("Number of generations")
    ylabel("Frequency of altruists")

    # Save fig
    savefig("freqA_b=$(b)_beta=$beta.png")

    # Close fig
    close()

end

#PRINT RESULTS TO A FILE
function Save_file(x, y, filename)
    arq = open(filename,"w")

    for cont in 1:length(x)
        write(arq,string(x[cont]))
        write(arq,"\t")
        write(arq,string(y[cont]))
        write("\n")
    end

    close(arq)
end


# Funcao rho: calcula (autovalor dominante - 1) a partir dos parametros
function rho(m,w0,b,c,n,deltaf, deltac ,beta)

    # Cria matriz de selecao
    S = matriz_selecao(n,deltaf, deltac ,b,c,w0,beta)

    # Cria matriz de migracao
    M = matriz_migracao(n,m)

    Prod = np.dot(M,S)
            
    return av_dominante(Prod)-1
end

# Encontra m critico ate o qual todos os autovalores sao maiores que 1 (se 0<m<1, usa o metodo da secante para achar raizes)
function m_critico(w0,b,c,n,deltaf, deltac,beta)

    fa = rho(0.,w0,b,c,n,deltaf, deltac,beta)

    fb = rho(1.,w0,b,c,n,deltaf, deltac ,beta)

    if(fa*fb>0)
        if(fa>0)
            return 1.
        else
            return 0.
        end
    end

    # brenth encontra o zero da funcao rho no intervalo [0.,1.]
    return so.brentq(rho, 0., 1., xtol=0.001, args=(w0,b,c,n,deltaf, deltac,beta))
end

# Fitness - recebe o tipo (0 ou 1 - altruista ou nao altruista)
# o numero de altruistas no grupo e o numero total de individuos no grupo
function fitness(tipo,k,n,delta,b,c,w0)
    return w0 + delta*(b*(k/n) - (1-tipo)*c)
end
    
# Fitness medio - recebe k e n
function fitness_m(k,n,delta,b,c,w0)
    freq_A = k/n

    return freq_A*fitness(0,k,n,delta,b,c,w0)+(1-freq_A)*fitness(1,k,n,delta,b,c,w0)
end

# Probabilidade de vencer - recebe k e n
# calcula a probabilidade de vencer um combate
function pVencer(deltac,k,n)
    return 1/(1+exp(-4*deltac*k/n))
end
    
# Select - recebe indice k,l e o numero total n
# calcula o numero medio de grupos do tipo k que viram tipo l por força de selecao
function sel(k,l,n,deltaf,deltac,b,c,w0,beta)
    w = fitness(0,k,n,deltaf,b,c,w0)
    p = k*w/(n*fitness_m(k,n,deltaf,b,c,w0))
    
    d = Binomial(n,p)
    
    return (1-beta+2*beta*pVencer(deltac,k,n))*pdf(d,l)
end
    
    
# Matriz de selecao - recebe n e parametros delta,b,c,w0 e beta
function matriz_selecao(n,deltaf,deltac,b,c,w0,beta)

    #S = np.fromfunction((i,j)-> sel(i,j,n,deltaf,deltac,b,c,w0,beta),(n+1,n+1))
    
    S = Array(Float64, n+1,n+1)
    
    for k in 1:n+1
        for l in 1:n+1
            S[(n+1)*(l-1)+k] = sel(k-1,l-1,n,deltaf,deltac,b,c,w0,beta)
        end
    end
    
    return S
end    
    
# Mig - recebe indice k,l e o numero total n
# Calcula probabilidade de um grupo k virar tipo l em relacao ao efeito de migracao
function mig(k,l,m)

    d = Binomial(k,1-m)

    #return sm.comb(k,l)*((1-m)^l)*(m^(k-l))
    return pdf(d,l)
end
     
# Matriz do numero medio de grupos que eram do tipo k e viram tipo l depois da migracao
function matriz_mig1(n,m)

    #M = np.fromfunction((i,j)-> mig(i,j,m),(n+1,n+1))
    M = Array(Float64, n+1, n+1)
    
    for k in 1:n+1
        for l in 1:n+1
            M[(n+1)*(l-1)+k] = mig(k-1,l-1,m)
        end
    end
    
    return M
end    
    
# Calcula o numero medio de grupos que recebem um altruista (eram do tipo 0 e viram do tipo 1)
function func(k,m)

    return m*k
end
    
# Matriz do numero medio de grupos que recebem altruistas
function matriz_mig2(n,m)

    T = zeros(n+1,n+1)

    T[:,2] = func(np.arange(n+1),m)

    return T
end
    
# Matriz de migracao - recebe n
function matriz_migracao(n,m)

    M = (m==0)? np.identity(n+1) : matriz_mig1(n,m)

    T = matriz_mig2(n,m)

    return M+T
end

# Calcula autovalor dominante
function av_dominante(Matriz)

    Matriz = np.delete(Matriz,0,0)

    Matriz = np.delete(Matriz,0,1)

    av = sl.eigvals(Matriz)
    
    x = sortperm(av, by=abs, rev=true)
    
    ind = x[1]
    
    return np.real(av[ind])[1]
end
