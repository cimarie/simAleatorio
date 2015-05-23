using Distributions
using PyCall
@pyimport smtplib
@pyimport collections
@pyimport itertools

# Main
function main()

    n = 26 #NUMBER OF MEMBERS IN EACH DEME
    N = 5000 #NUMBER OF DEMES
    A = 3  #INITIAL NUMBER OF ALTRUISTS
    beta = 0.8 #THE PROBABILITY OF A DEME TO ENTER IN A CONTEST
    LA = 1. #PARAMETER THAT INFLUENCES THE PROBABILITY OF WINNING (L)
#    b = 0.15 #BENEFIT FROM AN INTERACTION
    b = 0.
    c = 10. #COST TO THE ALTRUIST IN THE INTERACTION
    delta = 0.01 #SELECTION STRENGTH
    PM = 0.0001 #MUTATION PROBABILITY
    PMIG = 0.1 #MIGRATION PROBABILITY

    #vec_d = [.01, .1, .5, .8, 1.]
    vec_d = [delta]

    texto = "Simulacao: b=$b, c=$c, lambda=$LA, beta=$beta, mu=$PM, N=$N, n=$n\n\n"

    for delta in vec_d

        println("O tempo ate o fim eh")
        it, ultimo, itmeio, meio = Simulate(N,n,A,LA,b,c,beta,delta,PM,PMIG)

        texto =  texto * "Para delta = $delta, o tempo ate o fim eh: $(it+1) e a frequencia de altruistas inicial eh $(N*n/(N*n)) e a frequencia final de altruistas eh $ultimo\n  (quando t = $itmeio, freq = $meio\n\n)"

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

    IT = 7000
    precisao = 0.01
    stop = 0

    crit = A == N*n ? 0. : 1.

    #set counter zero
    it = 0

    if A < N*n
        #GENERATES THE DEMES MATRIX WITH THE GIVEN PARAMETERS OF N and n:
        DEMES = zeros(Int32,N,n)

        #println(DEMES)

        DEMES = Initialize(N,n,A,DEMES)
    else 
        DEMES = ones(Int32,N,n)
    end

    #println(DEMES)

    #ALTRUISTS COUNTER:
    AL = Float64[]
    push!(AL,(sum(DEMES))/(N*n))

    while it < IT
    
        #CONDITION TO STOP:
        if abs(AL[it+1]-crit)<precisao
			println("Acabou na geracao $it")
            break
        end
        
        DEMES = Make_contest(DEMES, N, n, LA, b, c, beta, delta)
        DEMES = Insert_mutation(N, n, PM, DEMES)
        DEMES = Make_migration(N, n, PMIG, DEMES)

        #ALTRUISTS COUNTER:
        push!(AL,(sum(DEMES))/(N*n))
        println(sum(DEMES)/(N*n))

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
    
    return it-1, AL[it], it/2, AL[it/2]
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

    if beta == 0

        #Group reproduction
        knums = [sum(line) for line in DEMES]
        distk = collections.Counter(knums)
        xk = keys(distk)
        pk = [i*lfitm[i] for i in xk]
        pk = pk/sum(pk)

        #
    else

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
            L = (exp(4*LA*delta*xD))/(exp(4*LA*delta*xD)+exp(4*LA*delta*xd))
            ##CALCULATE PA, THE PROBABILITY OF AN INDIVIDUAL (D or d) BE AN ALTRUIST:
            W_AD = 1 + delta*(-c + b*xD) #FITNESS OF D ALTRUISTS
            W_Ad = 1 + delta*(-c + b*xd) #FITNESS OF d ALTRUISTS
            W_ND = 1 + delta*b*xD #FITNESS OF D NON-ALTRUISTS
            W_Nd = 1 + delta*b*xd #FITNESS OF d NON-ALTRUISTS
            W_mD = W_AD*xD + W_ND*(1-xD) #AVERAGE FITNESS OF D
            W_md = W_Ad*xd + W_Nd*(1-xd) #AVERAGE FITNESS OF d

            PAD = xD*W_AD/W_mD
            PAd = xd*W_Ad/W_md

            taxaD = PAD
            taxad = PAd

            ## THE CONTEST HAPPENS: WINNER TAKES 2 DEMES, LOSER TAKES 0 DEMES
            if rand() < beta
            
                if rand() < L
                    taxad = PAD
                else
                    taxaD = PAd            
                end
                
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

    end

    return DEMES
end

# Conflito entre os grupos
function Conflito(N, n, knums, beta, lfit, lfitm, mpvencer)

    if beta == 0
        return Reproducao_grupo(N, n, knums, lfit, lfitm)

    else
        # Calcula numero de grupos que se envolvem em conflitos
        lconts = numpy.random.binomial(N,beta)
        indices = numpy.random.permutation(N)
        lconts = lconts if lconts%2==0 else lconts-1
        lista_conflitos = izip(indices[0:lconts/2], indices[lconts/2:lconts])

        for elem in lista_conflitos:
            i,j = elem
            ki,kj = knums[i], knums[j]
            venc,perd = (i,j) if mpvencer[ki,kj] > random.random() else (j,i)
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

function Make_migration(N, n, PMIG, DEMES)
    
    if PMIG == 0.
        return DEMES
    end
    
    #MIGRATION
    VDEME = Int32[]
    VPOS = Int32[]
    VMIG = Int32[]
    ##CHOOSES RANDOMLY SOME INDIVIDUALS FROM EACH DEME TO MIGRATE AND PUTS IN VMIG:
    for group in 1:N
        
        migrators = rand(Binomial(n, PMIG))
        
        if migrators > 0
            VPOS = vcat(VPOS, [1:migrators])
        end

        while migrators > 0
            # Take the first ones to migrate
            push!(VMIG, DEMES[group,migrators])
            # Fill VDEME with the groups that have blank spaces
            push!(VDEME,group)
            migrators -= 1
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
        write(arq, "\n")
    end

    close(arq)
end
