using Distributions
using PyCall
using Iterators
@pyimport smtplib
@pyimport collections
@pyimport itertools
@pyimport numpy as np
@pyimport scipy.stats as ss

# Main
function main()

    #n = 26 #NUMBER OF MEMBERS IN EACH DEME
    N, n = 1, 5000
    #N = 5000 #NUMBER OF DEMES
    A = N*n-1  #INITIAL NUMBER OF ALTRUISTS
    beta = 0. #THE PROBABILITY OF A DEME TO ENTER IN A CONTEST
    LA = 2. #PARAMETER THAT INFLUENCES THE PROBABILITY OF WINNING (L)
    b = 0. #BENEFIT
    c = 10. #COST TO THE ALTRUIST IN THE INTERACTION
    delta = 0.01 #SELECTION STRENGTH
    #PM = 0.0001 #MUTATION PROBABILITY
    PM = 0.
    PMIG = 0. #MIGRATION PROBABILITY

    #vec_d = [.01, .1, .5, .8, 1.]
    vec_d = [delta]

    texto = "Simulacao: b=$b, c=$c, lambda=$LA, beta=$beta, mu=$PM, N=$N, n=$n\n\n"

    for delta in vec_d

        println("O tempo ate o fim eh")
        Simulate(N,n,A,LA,b,c,beta,delta,PM,PMIG)

        #texto =  texto * "Para delta = $delta, o tempo ate o fim eh: $(it+1) e a frequencia de altruistas inicial eh $(N*n/(N*n)) e a frequencia final de altruistas eh $ultimo\n  (quando t = $itmeio, freq = $meio\n\n)"

    end

    #manda_email(texto)

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
function Simulate(N,n,A,LA,b,c,beta,delta,PM,PMIG)

    IT = 5000
    precisao = 0.01

    #set counter one
    it = 1

    DEMES = Initialize(N,n,A)

    #ALTRUISTS COUNTER:
    AL = Float64[]
    push!(AL,countnz(DEMES)/(N*n))

    lfit = cria_listaf(n, b, c, delta)
    lfitm = cria_listafm(n, b, c, delta)
    mpvencer = matriz_vencedores(LA, n)

    crit = AL[1] == 1. ? 0. : 1.
    
    while it < IT
    
        #CONDITION TO STOP:
        if abs(AL[it]-crit)<precisao
			println("Acabou na geracao $it")
            break
        end
        
        knums = [countnz(DEMES[i,:]) for i in [1:N]]
        glabels = Make_contest(N, knums, beta, lfitm, mpvencer)
        DEMES = Indiv_reprod(N, n, lfit, lfitm, glabels)
        DEMES = Insert_mutation(N, n, PM, DEMES)
        DEMES = Make_migration(N, n, PMIG, DEMES)

        #ALTRUISTS COUNTER:
        push!(AL,countnz(DEMES)/(N*n))
        #println(countnz(DEMES)/(N*n))
        @printf "Na geracao %d ------> %.7f\n " it countnz(DEMES)/(N*n)

        it += 1
    end
    
    #ARRAY WITH THE PERIODS OF TIME FOR WHICH # ALTRUISMS WAS MEASURED
    IC = [0:it]

    println(AL[1])
    println(AL[it])

    #SAVE RESULTS TO A FILE
    #Save_file(IC, AL,"sim.txt")

    # Save results to a figure
    #Save_fig(N, n, A, IC, AL, b, beta, PMIG)
    
    #return it-1, AL[it], it/2, AL[it/2]
end

function cria_listaf(n, b, c, delta)
    return [fitness(0,n,k,b,c,delta) for k in [0:n]]
end

function cria_listafm(n, b, c, delta)
    return [fitness_m(n,k,b,c,delta) for k in [0:n]]
end

function fitness(tipo, n, k, b, c, delta)
    return 1 + delta*(b*k/n - (1-tipo)*c)
end

function fitness_m(n, k, b, c, delta)
    freq = k/n
    return freq*fitness(0,n,k,b,c,delta)+(1-freq)*fitness(1,n,k,b,c,delta)
end

# Probabilidade de i vencer j
function pVencer(i, j, alpha, n)
    arg1 = 4.*alpha*i/n
    arg2 = 4.*alpha*j/n
    return exp(arg1)/(exp(arg1)+exp(arg2))
end

# Constroi matriz com as probabilidades de vencer (linha vencer coluna)
function matriz_vencedores(alpha, n)
    w = [0:n]
    v = reshape(w,1,n+1)
    c = np.meshgrid(v,w)
    return np.map((i,j)->pVencer(i,j,alpha,n),c[2],c[1])
    #return broadcast((i,j)->pVencer(i,j,alpha,n), w, v)
end

function Regression(N,n,A,LA,b,c,beta,delta,PM)

    it = Int32[] 

	vec_m = [0:0.1:1]

	for PMIG in vec_m
	    #println(PMIG)
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
function Initialize(n_groups, n_ind, n_A)

    if n_A == n_groups*n_ind
        return ones(Int32, (n_groups,n_ind))
    end

    DEMES = Array(Float64, n_groups, n_ind)
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

function Group_reprod(N, knums, lfitm)

    #Group reproduction
    distk = collections.Counter(knums)
    xk = keys(distk)
    xk = Int32[x for x in xk]
    pk = Float64[distk[i]*lfitm[i+1] for i in xk]
    pk = pk/sum(pk)
    distw = Categorical(pk)
    glabels = rand(distw, N)
    glabels = Int32[xk[i] for i in glabels]

    return glabels
end

# Reproducao individual (ocorre apos o conflito/reproducao de grupo)
function Indiv_reprod(N, n, lfit, lfitm, glabels)

    # Reproducao individual
    aux = zeros(N,n)

    for i in [1:N]
        k = glabels[i] 
        wA = lfit[k+1] 
        wm = lfitm[k+1]
        filhosA = rand(Binomial(n, k*wA/(n*wm)))
        aux[i,1:filhosA] = 1
        aux[i,:] = shuffle(vec(aux[i,:]))
    end

    return aux 
end

# Make the contest between two demes
function Make_contest(N, knums, beta, lfitm, mpvencer)

    if N == 1
        return knums 
    end

    if beta == 0
        return Group_reprod(N,knums,lfitm)

    else
        # Calcula numero de grupos que se envolvem em conflitos
        lconts = rand(Binomial(N,beta))
        indices = shuffle([1:N])
        lconts = mod(lconts,2)==0? lconts : lconts-1
        lista_conflitos = itertools.izip(indices[1:lconts/2], indices[lconts/2+1:lconts])

        for elem in lista_conflitos
            i,j = elem
            ki,kj = knums[i], knums[j]
            venc,perd = mpvencer[ki+1,kj+1] > rand() ? (i,j) : (j,i)
            knums[perd]=knums[venc]
        end
    end

    return knums
end

function Insert_mutation(n_groups, n_ind, mu, DEMES)

    #MUTATION
    flips = rand!(Bernoulli(mu),Array(Int32,n_groups,n_ind))
    DEMES = mod(DEMES+flips,2)

    return DEMES
end

# Migracao entre os grupos
function Make_migration(N, n, pmig, DEMES)
    
    if pmig > 0.

        # Lista com numero de migrantes por grupo
        nmigs = rand(Binomial(n,pmig),N)

        # Lista com os migrantes
        migs = Int32[]
        for i in 1:N
            for j in 1:nmigs[i]
                push!(migs, DEMES[i,j])
            end
        end
        migs = shuffle(migs)

        cont = 1
        for i in 1:N
            DEMES[i, 1:nmigs[i]] = migs[cont:cont+nmigs[i]-1]
            cont += nmigs[i]
            DEMES[i,:] = shuffle(vec(DEMES[i,:]))
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
        write(arq, "\n")
    end

    close(arq)
end
