using DataFrames
using LinearAlgebra
using Optim

function threepoint(tree::AbstractTree, trait::String, N::Int64, nodes)
    #prefrom algortihm in Ho & Ane 2014
    #function estimaterates gets inputs into right form
    #nodes - vector of nodes in the traversal order

    for i in 1:N
        #need to see if node is a tip (leaf)
        if isleaf(tree, nodes[i])
    
            #calculations
            nodet = getnodedata(tree, nodes[i])["t"]
            nodetrait = getnodedata(tree, nodes[i])[trait]
            
            #update nodes
            setnodedata!(tree, nodes[i], "logV", log(nodet))
            setnodedata!(tree, nodes[i], "p", 1/nodet)
            setnodedata!(tree, nodes[i], "yl", nodetrait)
            setnodedata!(tree, nodes[i], "xl", 1)
            setnodedata!(tree, nodes[i], "Q", nodetrait / nodet)
            setnodedata!(tree, nodes[i], "xx", 1/nodet)
            setnodedata!(tree, nodes[i], "yy", nodetrait * nodetrait / nodet)
    
        else
    
            #need to find direct desendents 
            children = getchildren(tree, nodes[i])
    
            #child data
            childdata = getnodedata.(tree, children)
    
            #get needs the name of the vector of dicts, what we're looking for in each dict, and what to say if it cant find anything
            childp = get.(childdata, "p", missing) 
            childlogV = get.(childdata, "logV", missing)
            childxl = get.(childdata, "xl", missing)
            childyl = get.(childdata, "yl", missing)
            childQ = get.(childdata, "Q", missing)
            childxx = get.(childdata, "xx", missing)
            childyy = get.(childdata, "yy", missing)
    
            #nodedata
            nodet = getnodedata(tree, nodes[i])["t"]
     
            #calculations
            pA = sum(childp)
            ws = childp / pA
           
            newlogV = sum(childlogV) + log(1.0 + nodet * pA)
            newp = pA / (1.0 + nodet * pA)
            newxl = ws ⋅ childxl
            newyl = ws ⋅ childyl
            newQ = sum(childQ) - ((nodet * pA^2) / (1.0 + nodet * pA)) * newxl * newyl
            newxx = sum(childxx) - ((nodet * pA^2) / (1.0 + nodet * pA)) * newxl * newxl
            newyy = sum(childyy) - ((nodet * pA^2) / (1.0 + nodet * pA)) * newyl * newyl 
    
            #update nodes
            setnodedata!(tree, nodes[i], "logV", newlogV)
            setnodedata!(tree, nodes[i], "p", newp)
            setnodedata!(tree, nodes[i], "yl", newyl)
            setnodedata!(tree, nodes[i], "xl", newxl)
            setnodedata!(tree, nodes[i], "Q", newQ)
            setnodedata!(tree, nodes[i], "xx", newxx)
            setnodedata!(tree, nodes[i], "yy", newyy)
    
        end
    end

    return tree
end

function estimaterates!(tree::AbstractTree, trait::String)

    #get information from tree in order to preform threepoint
    nodes = getnodes(tree, postorder)
    N = length(nodes)
    n = length(getleaves(tree))

    tree = threepoint(tree, trait, N, nodes)

    #information from last node
    nodeN = getnodedata(tree, nodes[N])

    betahat = inv(nodeN["xx"]) * nodeN["Q"]
    sigmahat = (nodeN["yy"] - 2 * betahat * nodeN["Q"] + betahat * nodeN["xx"] * betahat) / n


    if sigmahat < 0 
        leaves = getleaves(tree)
        leafdata = getnodedata.(tree, leaves)
        leavestraits = get.(leafdata, "trait", missing)

        newtraits = ones(n)*betahat - leavestraits
        setnodedata!.(tree, leaves, "trait", newtraits)
        tree2 = threepoint(tree, trait, N, nodes)
        sigmahat = getnodedata(tree2, nodes[N])["yy"] / n
    end


    
    negloglik = (1/2) * (n * log(2π) + nodeN["logV"] + n + n * log(sigmahat))

    return betahat, sigmahat, negloglik, tree
end

function estimaterates(tree::AbstractTree, trait::String)
    #Returns evolution rate, starting value and negative log loglikelihood for traits on tip of tree
    #INPUTS
    #tree = tree with lengths, leaves all same length, trait data on leaves
    #trait = string with name of trait as found on leaves
    #OUTPUTS
    #sigmahat - evolution rate
    #betahat - estimated root trait value
    #negloglik - negative loglikelihood
    #tree - tree w/ intermediate calculations on nodes

    #need to add meaningful error when cant find traits on tree leaves

    nodes = getnodes(tree, postorder)
    N = length(nodes)

    setnodedata!.(tree, nodes, "t", undef) #could do this as I create it I guess
    setnodedata!.(tree, nodes, "logV", undef)
    setnodedata!.(tree, nodes, "p", undef)
    setnodedata!.(tree, nodes, "Q", undef)
    setnodedata!.(tree, nodes, "xl", undef)
    setnodedata!.(tree, nodes, "xx", undef)
    setnodedata!.(tree, nodes, "yl", undef)
    setnodedata!.(tree, nodes, "yy", undef)

    for i in 1:(N-1)
        par = getinbound(tree, nodes[i])
        length = getlength(tree, par)
        setnodedata!(tree, nodes[i], "t", length)
    end

    setnodedata!(tree, nodes[N], "t", 0)

    return estimaterates!(tree, trait)

end

######################################################################################################
#3 Point Signal

function estimaterates(tree::AbstractTree, trait::String, lambda::Float64)
    #Returns evolution rate, starting value and negative log loglikelihood for traits on tip of tree
    #INPUTS
    #tree = tree with lengths, leaves all same length, trait data on leaves
    #trait = string with name of trait as found on leaves
    #OUTPUTS
    #sigmahat - evolution rate
    #betahat - estimated root trait value
    #negloglik - negative loglikelihoods
    #tree - tree w/ intermediate calculations on nodes

    #need to add meaningful error when cant find traits on tree leaves

    nodes = getnodes(tree, postorder)
    N = length(nodes)

    setnodedata!.(tree, nodes, "t", undef) #could do this as I create it I guess
    setnodedata!.(tree, nodes, "logV", undef)
    setnodedata!.(tree, nodes, "p", undef)
    setnodedata!.(tree, nodes, "Q", undef)
    setnodedata!.(tree, nodes, "xl", undef)
    setnodedata!.(tree, nodes, "xx", undef)
    setnodedata!.(tree, nodes, "yl", undef)
    setnodedata!.(tree, nodes, "yy", undef)


    for i in 1:(N-1)
        par = getinbound(tree, nodes[i])
        length = getlength(tree, par)
        setnodedata!(tree, nodes[i], "t", length)
    end

    setnodedata!(tree, nodes[N], "t", 0)

    t = getnodedata.(tree, nodes, "t")

    #=
    #for non leaves mult by lambda - could do this w/o a loop if get list of internal nodes
    for i in 1:N
        if !isleaf(tree, nodes[i])
            t = lambda * getnodedata(tree, nodes[i], "t")
            setnodedata!(tree, nodes[i], "t", t)
        end
    end
    =#

    return estimaterates!(tree, trait, lambda, t)

end

function tooptimise(lambda, tree, nodes, t, trait, N, n)
    #lambda - value for Signal   
    #tree - tree
    #nodes - nodes of the tree in traversal order
    #t - vector of og branch lengths
    #trait - string of what trait is called on the tree
    #N - total number of nodes
    #n - number of leaves

    #for non leaves mult by lambda - could do this w/o a loop if get list of internal nodes
    for i in 1:N
        if !isleaf(tree, nodes[i])
            tupdate = lambda[1] * t[i]
            setnodedata!(tree, nodes[i], "t", tupdate) #tupdate is being put through as a vector
        end
    end

    tree = threepoint(tree, trait, N, nodes)

    #information from last node
    nodeN = getnodedata(tree, nodes[N])

    betahat = inv(nodeN["xx"]) * nodeN["Q"]
    sigmahat = (nodeN["yy"] - 2 * betahat * nodeN["Q"] + betahat * nodeN["xx"] * betahat) / n
    
    
    if sigmahat < 0 
        leaves = getleaves(tree)
        leafdata = getnodedata.(tree, leaves)
        leavestraits = get.(leafdata, "trait", missing)
    
        newtraits = ones(n)*betahat - leavestraits
        setnodedata!.(tree, leaves, "trait", newtraits)
        tree2 = threepoint(tree, trait, N, nodes)
        sigmahat = getnodedata(tree2, nodes[N])["yy"] / n
    end
    
    
        
    negloglik = (1/2) * (n * log(2π) + nodeN["logV"] + n + n * log(sigmahat))

    return negloglik
end

function estimaterates!(tree::AbstractTree, trait::String, lambda::Float64, t)

    nodes = getnodes(tree, postorder)
    N = length(nodes)
    n = length(getleaves(tree))

    #info for optimiser
    lower = [0.0]
    upper = [1.0]
    start = [lambda]

    #optimise to find lambda
    opts = optimize(x -> tooptimise(x, tree, nodes, t, trait, N, n), lower, upper, start, Fminbox(LBFGS())) 
    lambda = Optim.minimizer(opts)

    #update internal branches
    for i in 1:N
        if !isleaf(tree, nodes[i])
            tupdate = lambda[1] * t[i]
            setnodedata!(tree, nodes[i], "t", tupdate) #tupdate is being put through as a vector
        end
    end

    #preform threepoint on lambda transformed tree
    tree = threepoint(tree, trait, N, nodes)

    #information from last node
    nodeN = getnodedata(tree, nodes[N])

    betahat = inv(nodeN["xx"]) * nodeN["Q"]
    sigmahat = (nodeN["yy"] - 2 * betahat * nodeN["Q"] + betahat * nodeN["xx"] * betahat) / n


    if sigmahat < 0 
        leaves = getleaves(tree)
        leafdata = getnodedata.(tree, leaves)
        leavestraits = get.(leafdata, "trait", missing)

        newtraits = ones(n)*betahat - leavestraits
        setnodedata!.(tree, leaves, "trait", newtraits)
        tree2 = threepoint(tree, trait, N, nodes)
        sigmahat = getnodedata(tree2, nodes[N])["yy"] / n
    end

    negloglik = (1/2) * (n * log(2π) + nodeN["logV"] + n + n * log(sigmahat))

    return lambda[1], betahat, sigmahat, negloglik, tree
end

#####################################################################################################################
#=
function threepoint(tree::AbstractTree, df::DataFrame, traits::DataFrame)
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

    df = threepoint(tree, df, traits)    

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

    return estimaterates!(tree, traits, df)
end
=#
