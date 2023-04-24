using DataFrames
using LinearAlgebra
using Optim

mutable struct TraitData
    name::String
    value::Float64
    t::Float64
    logV::Float64
    p::Float64
    yl::Float64
    xl::Float64
    Q::Float64
    xx::Float64
    yy::Float64
end

traitdata(name::String, value::Float64, t::Float64 = 0.0) =
    TraitData(name, value, t, NaN, NaN, NaN, NaN, NaN, NaN, NaN)

TraitData() = traitdata("", NaN)

const TLB{RT, LenUnits} = LinkBranch{RT, String, Nothing, LenUnits}
const TLN{RT, LenUnits} = LinkNode{RT, String, TraitData, TLB{RT, LenUnits}}
const TLT{RT, TD, LenUnits} = LinkTree{RT, String, TLN{RT, LenUnits}, TLB{RT, LenUnits}, TD}
const TLTD{RT, LenUnits} = TLT{RT, Dict{String, Any}, LenUnits}
const TraitTree = TLTD{OneRoot, Float64}

function threepoint!(tree::T, _::String, nodes::Vector{N}) where 
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}}
    #prefrom algortihm in Ho & Ane 2014
    #function estimaterates gets inputs into right form
    #nodes - vector of nodes in the traversal order

    for node in nodes
        nd = getnodedata(tree, node)
        nodet = nd.t

        #need to see if node is a tip (leaf)
        if isleaf(tree, node)
            nodetrait = nd.value
            
            # update node data
            nd.logV = log(nodet)
            nd.p = inv(nodet)
            nd.yl = nodetrait
            nd.xl = 1.0
            nd.Q = nodetrait / nodet
            nd.xx = inv(nodet)
            nd.yy = nodetrait * nodetrait / nodet
        else
            # need to find direct desendents 
            children = getchildren(tree, node)

            # child data
            childdata = getnodedata.(tree, children)
    
            # calculations
            childp = [s.p for s in childdata]
            pA = sum(childp)
            ws = childp / pA
           
            calc = 1.0 + nodet * pA
            nd.logV = sum(s.logV for s in childdata) + log(calc)
            nd.p = pA / calc
            nd.xl = ws ⋅ [s.xl for s in childdata]
            nd.yl = ws ⋅ [s.yl for s in childdata]

            c2 = nodet * pA^2 / calc
            nd.Q = sum(s.Q for s in childdata) - c2 * nd.xl * nd.yl
            nd.xx = sum(s.xx for s in childdata) - c2 * nd.xl * nd.xl
            nd.yy = sum(s.yy for s in childdata) - c2 * nd.yl * nd.yl
        end
        # @assert getnodedata(tree, node) === nd
    end
end

function estimaterates!(tree::T, trait::String) where T <: AbstractTree

    #get information from tree in order to preform threepoint
    nodes = getnodes(tree, postorder)
    n = nleaves(tree)

    threepoint!(tree, trait, nodes)

    #information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q
    sigmahat = (nd.yy - 2 * betahat * nd.Q + betahat * nd.xx * betahat) / n

    if sigmahat < 0 
        leaves = getleaves(tree)
        for leaf in leaves
            ld = getnodedata(tree, leaf)
            ld.value = betahat - ld.value
        end

        threepoint!(tree, trait, nodes)
        sigmahat = nd.yy / n
    end
    
    negloglik = (1.0 / 2.0) * (n * log(2π) + nd.logV + n + n * log(sigmahat))

    return betahat, sigmahat, negloglik
end

function estimaterates(tree::T, trait::String) where T <: AbstractTree
    #Returns evolution rate, starting value and negative log loglikelihood for traits on tip of tree
    #INPUTS
    #tree = tree with lengths, leaves all same length, trait data on leaves
    #trait = string with name of trait as found on leaves
    #OUTPUTS
    #sigmahat - evolution rate
    #betahat - estimated root trait value
    #negloglik - negative loglikelihood

    #need to add meaningful error when cant find traits on tree leaves

    nodes = getnodes(tree, anyorder)

    for node in nodes
        val = getnodedata(tree, node).value
        if hasinbound(tree, node)
            len = _getlength(tree, _getinbound(tree, node))
            td = traitdata(trait, val, len)
            setnodedata!(tree, node, td)
        else
            td = traitdata(trait, val)
            setnodedata!(tree, node, td)
        end
    end

    return estimaterates!(tree, trait)
end

######################################################################################################
#3 Point Signal

function estimaterates(tree::T, trait::String, lambda::Float64) where T <: AbstractTree
    #Returns evolution rate, starting value and negative log loglikelihood for traits on tip of tree
    #INPUTS
    #tree = tree with lengths, leaves all same length, trait data on leaves
    #trait = string with name of trait as found on leaves
    #OUTPUTS
    #sigmahat - evolution rate
    #betahat - estimated root trait value
    #negloglik - negative loglikelihoods

    #need to add meaningful error when cant find traits on tree leaves

    nodes = getnodes(tree, postorder)

    for node in nodes
        val = getnodedata(tree, node).value
        if hasinbound(tree, node)
            len = _getlength(tree, _getinbound(tree, node))
            setnodedata!(tree, node, traitdata(trait, val, len))
        else
            setnodedata!(tree, node,
                         traitdata(trait, getnodedata(tree, node).value))
        end
    end

    t = [getnodedata(tree, node).t for node in nodes]

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

function tooptimise(lambda::Vector{Float64}, tree::T, nodes::Vector{N},
                    t::Vector{Float64}, trait::String, n::Int) where 
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}}
    #lambda - value for Signal   
    #tree - tree
    #nodes - nodes of the tree in traversal order
    #t - vector of og branch lengths
    #trait - string of what trait is called on the tree
    #N - total number of nodes
    #n - number of leaves

    #for non leaves mult by lambda - could do this w/o a loop if get list of internal nodes
    for (i, node) in enumerate(nodes)
        if !isleaf(tree, node)
            nd = getnodedata(tree, node)
            nd.t = lambda[1] * t[i]
        end
    end

    threepoint!(tree, trait, nodes)

    #information from last node
    nd = getnodedata(tree, last(nodes))

    betahat = inv(nd.xx) * nd.Q
    sigmahat = (nd.yy - 2.0 * betahat * nd.Q + betahat * nd.xx * betahat) / n  
    
    if sigmahat < 0.0 
        leaves = getleaves(tree)
        for leaf in leaves
            nd = getnodedata(tree, leaf)
            nd.value = betahat - nd.value
        end

        threepoint!(tree, trait, nodes)
        sigmahat = nd.yy / n
    end

    negloglik = (1/2) * (n * log(2π) + nd.logV + n + n * log(sigmahat))

    return negloglik
end

function estimaterates!(tree::T, trait::String, lambda::Float64, t::Vector{Float64}) where T <: AbstractTree

    nodes = getnodes(tree, postorder)
    N = length(nodes)
    n = length(getleaves(tree))
    root = nodes[N]

    #info for optimiser
    lower = [0.0]
    upper = [1.0]
    start = [lambda]

    #optimise to find lambda
    opts = optimize(x -> tooptimise(x, tree, nodes, t, trait, n),
                    lower, upper, start, Fminbox(LBFGS()))
    lambda = Optim.minimizer(opts)

    #update internal branches
    for (i, node) in enumerate(nodes)
        if !isleaf(tree, node)
            tupdate = lambda[1] * t[i]
            getnodedata(tree, node).t = tupdate # tupdate is being put through as a vector
        end
    end

    #preform threepoint on lambda transformed tree
    threepoint!(tree, trait, nodes)

    #information from last node
    nd = getnodedata(tree, root)

    betahat = inv(nd.xx) * nd.Q
    sigmahat = (nd.yy - 2 * betahat * nd.Q + betahat * nd.xx * betahat) / n

    if sigmahat < 0 
        leaves = getleaves(tree)
        leafdata = getnodedata.(tree, leaves)
        leavestraits = get.(leafdata, "trait", missing)

        newtraits = betahat .- leavestraits
        getnodedata.(tree, leaves, "trait", newtraits)
        threepoint!(tree, trait, nodes)
        sigmahat = getnodedata(tree, root).yy / n
    end

    negloglik = (1.0 / 2.0) * (n * log(2π) + nd.logV + n + n * log(sigmahat))

    return lambda[1], betahat, sigmahat, negloglik
end
