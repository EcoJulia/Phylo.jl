using DataFrames

function threepoint(tree, df, traits, N)
    #prefrom algortihm in Ho & Ane 2014
    #function estimaterates gets inputs into right form
    for i in 1:N
        #need to see if node is a tip (leaf)
        if isleaf(tree, df.Node[i])

            #where this leaf is in traits df and save trait
            data = filter(:species => x -> x == df.Node[i], traits)[1,2]
    
            df.logV[i] = log(df.t[i]) 
            df.p[i] = 1/df.t[i]
            df.yl[i] = data
            df.xl[i] = 1
            df.Q[i] = data/df.t[i]
            df.xx[i] = 1/df.t[i]
            df.yy[i] =  (data*data)/df.t[i]
        else
            #need to find direct desendents 
            children = getchildren(tree, df.Node[i])
            cdf = filter(:Node => x -> x ∈ children, df)
    
            pA = sum(cdf.p)
            ws = cdf.p/pA
    
            df.logV[i] = sum(cdf.logV) + log(1+ df.t[i]*pA) 
            df.p[i] = pA/(1+df.t[i]*pA)
            df.xl[i] = sum(ws.*cdf.xl)
            df.yl[i] = sum(ws.*cdf.yl)
            df.Q[i] = sum(cdf.Q) - ((df.t[i]*pA^2)/(1+df.t[i]*pA))*df.xl[i]*df.yl[i]
            df.xx[i] = sum(cdf.xx) - ((df.t[i]*pA^2)/(1+df.t[i]*pA))*df.xl[i]*df.xl[i]
            df.yy[i] = sum(cdf.yy) - ((df.t[i]*pA^2)/(1+df.t[i]*pA))*df.yl[i]*df.yl[i]
        end
    end
    return df
end


function estimaterates(tree, traits)
    #Returns evolution rate, starting value and negative log loglikelihood for traits on tip of tree
    #INPUTS
    #tree = tree with lengths, leaves all same length
    #traits = dataframe with leaf names and trait values, leafnames called 'species' and trait information called 'data'
    #OUTPUTS
    #sigmahat - evolution rate
    #betahat - estimated root trait value
    #negloglik - negative loglikelihood
    #df - dataframe of results - will add what all mean later

    #INPUT TESTS
    #traits for each leaf?

    nodenames = getnodenames(tree, breadthfirst)
    #total number of nodes
    N = length(nodenames)
    #number of leaves
    n = nrow(traits)

    #dataframe to save information in
    df = DataFrame(Node = nodenames, 
            t = Vector{Union{Nothing, Float64}}(nothing, N),
            logV = Vector{Union{Nothing, Float64}}(nothing, N),
            p = Vector{Union{Nothing, Float64}}(nothing, N),
            Q = Vector{Union{Nothing, Float64}}(nothing, N),
            xl = Vector{Union{Nothing, Float64}}(nothing, N),
            xx = Vector{Union{Nothing, Float64}}(nothing, N),
            yl = Vector{Union{Nothing, Float64}}(nothing, N),
            yy = Vector{Union{Nothing, Float64}}(nothing, N))

    for i in 1:(N-1)
        par = getparent(tree, df.Node[i])
        df.t[i] = distance(tree, df.Node[i], par)
    end

    df.t[N] = 0

    df = threepoint(tree, df, traits, N)    

    betahat = inv(df.xx[N]) * df.Q[N]
    sigmahat = (df.yy[N] - 2 * betahat * df.Q[N] + betahat * df.xx[N] * betahat)/n

    if sigmahat < 0 #if used prints df2 at the end?
        resdl = ones(n)*betahat - traits.data #replaces y which is traits
        traits2 = DataFrame(species = traits.species, data = resdl)
        df2 = threepoint(tree, df, traits2, N)
        sigmahat = df2.yy[N]/n
    end
    
    negloglik = (1/2)*(n*log(2*pi) + df.logV[N] + n + n*log(sigmahat))

    return betahat, sigmahat, negloglik, df
end  
