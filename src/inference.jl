using DataFrames
using LinearAlgebra

function threepoint(tree::AbstractTree, df::DataFrame, traits::DataFrame, N::Int64)
    #prefrom algortihm in Ho & Ane 2014
    #function estimaterates gets inputs into right form
    for row in eachrow(df)
        #need to see if node is a tip (leaf)
        if isleaf(tree, row.Node)

            #where this leaf is in traits df and save trait
            data = filter(:species => x -> x == row.Node, traits)[1,2]::Float64 #red in profile plot
    
            row.logV = log(row.t)
            row.p = 1.0 / row.t
            row.yl = data
            row.xl = 1.0
            row.Q = data / row.t
            row.xx = 1.0 / row.t
            row.yy =  (data * data) / row.t
        else
            #need to find direct desendents 
            children = getchildren(tree, row.Node)
            cdf = filter(:Node => x -> x ∈ children, df)
    
            pA = sum(cdf.p)
            ws = cdf.p / pA
    
            row.logV = sum(cdf.logV) + log(1.0 + row.t*pA) 
            row.p = pA/(1.0 + row.t * pA)
            row.xl = ws ⋅ cdf.xl #red in profile plot
            row.yl = ws ⋅ cdf.yl
            row.Q = sum(cdf.Q) - ((row.t * pA^2) / (1.0 + row.t * pA)) * row.xl * row.yl #red in profile plot
            row.xx = sum(cdf.xx) - ((row.t * pA^2) / (1.0 + row.t * pA)) * row.xl * row.xl
            row.yy = sum(cdf.yy) - ((row.t * pA^2) / (1.0 + row.t * pA)) * row.yl * row.yl
        end
    end
    return df
end

function estimaterates!(tree::AbstractTree, traits::DataFrame, df::DataFrame)
    #Returns evolution rate, starting value and negative log loglikelihood for traits on tip of tree
    #INPUTS
    #tree = tree with lengths, leaves all same length
    #traits = dataframe with leaf names and trait values, leafnames called 'species' and trait information called 'data'
    #df = dataframe to store threepoint values in
    #OUTPUTS
    #sigmahat - evolution rate
    #betahat - estimated root trait value
    #negloglik - negative loglikelihood
    #df - dataframe of results - will add what all mean later


    #total number of nodes
    N = nrow(df)
    #number of leaves
    n = nrow(traits)

    df = threepoint(tree, df, traits, N)    

    betahat = inv(df.xx[N]) * df.Q[N]
    sigmahat = (df.yy[N] - 2 * betahat * df.Q[N] + betahat * df.xx[N] * betahat) / n

    if sigmahat < 0 #if used prints df2 at the end?
        resdl = ones(n)*betahat - traits.data #replaces y which is traits
        traits2 = DataFrame(species = traits.species, data = resdl)
        df2 = threepoint(tree, df, traits2, N)
        sigmahat = df2.yy[N] / n
    end
    
    negloglik = (1/2) * (n * log(2π) + df.logV[N] + n + n * log(sigmahat))

    return betahat, sigmahat, negloglik, df #red in profile plot
end


function estimaterates(tree::AbstractTree, traits::DataFrame)
    nodes = getnodes(tree, postorder)
    #total number of nodes
    N = length(nodes)

    #dataframe to save information in
    df = DataFrame(Node = nodes, #empty dataframe outside the function
            t = Vector{Float64}(undef, N),
            logV = Vector{Float64}(undef, N),
            p = Vector{Float64}(undef, N),
            Q = Vector{Float64}(undef, N),
            xl = Vector{Float64}(undef, N),
            xx = Vector{Float64}(undef, N),
            yl = Vector{Float64}(undef, N),
            yy = Vector{Float64}(undef, N))

    for i in 1:(N-1)
        par = getinbound(tree, df.Node[i])
        df.t[i] = getlength(tree, par)
    end

    df.t[N] = 0

    for _ in 1:99
        estimaterates!(tree, traits, df)
    end
    return estimaterates!(tree, traits, df)
end
