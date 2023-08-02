using DataFrames
using LinearAlgebra
using Optim
using StaticArrays


mutable struct TraitData{NTraits}
    name::Vector{String}          #SizedVector{NTraits, String}
    value::Vector{Float64}        #SizedVector{NTraits, Float64}
    t::Float64
    logV::Float64
    p::Float64
    yl::Vector{Float64}           #SizedVector{NTraits, Float64}
    xl::Float64
    Q::Vector{Float64}            #SizedVector{NTraits, Float64}
    xx::Float64
    yy::Matrix{Float64}  #SizedMatrix{NTraits, NTraits, Float64}
end

traitdata(name::Vector{String}, value::Vector{Float64}, t::Float64 = 0.0) =
    TraitData{length(name)}(name, value, t, NaN, NaN, fill(NaN, length(name)), NaN, 
                            fill(NaN, length(name)), NaN, fill(NaN, length(name), length(name)))

TraitData{NTraits}() where NTraits = traitdata(fill("", NTraits), fill(NaN, NTraits))

const TLB{RT, LenUnits} = LinkBranch{RT, String, Nothing, LenUnits}
const TLN{NT, RT, LenUnits} = LinkNode{RT, String, TraitData{NT}, TLB{RT, LenUnits}}
const TLT{NT, RT, TD, LenUnits} = LinkTree{RT, String, TLN{NT, RT, LenUnits}, TLB{RT, LenUnits}, TD}
const TLTD{NT, RT, LenUnits} = TLT{NT, RT, Dict{String, Any}, LenUnits}
const TraitTree{NTraits} = TLTD{NTraits, OneRoot, Float64}


#=
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
=#

function threepoint!(tree::T, _::Vector{String}, nodes::Vector{N}) where 
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
            nd.yy = nodetrait * nodetrait' / nodet
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
            nd.xl = sum(ws .* [s.xl for s in childdata]) 
            nd.yl = sum(ws .* [s.yl for s in childdata])#, dims = 1) 

            c2 = nodet * pA^2 / calc
            nd.Q = sum(s.Q for s in childdata) - c2 * nd.xl * nd.yl
            nd.xx = sum(s.xx for s in childdata) - c2 * nd.xl * nd.xl
            nd.yy = sum(s.yy for s in childdata) - c2 * nd.yl * nd.yl'
        end
        # @assert getnodedata(tree, node) === nd
    end
end

function estimaterates!(tree::T, trait::Vector{String}, lambda) where T <: AbstractTree

    #get information from tree in order to preform threepoint
    nodes = getnodes(tree, postorder)
    n = nleaves(tree)

    if lambda !== missing
        t = [getnodedata(tree, node).t for node in nodes]

        #info for optimiser
        lower = [0.0]
        start = [lambda]

        #upper is going to be largest lambda can be without going past the leaves
        #want to find longest internal branch
        intnodeheights = nodeheights(tree, noleaves=true)
        longnodeheight = maximum(intnodeheights)

        leafnodeheights = nodeheights(tree, onlyleaves=true)
        shortleafheight = minimum(leafnodeheights)

        upper = [shortleafheight/longnodeheight]

        #optimise to find lambda
        opts = optimize(x -> tooptimise(x, tree, nodes, trait, n),
                        lower, upper, start, Fminbox(LBFGS()))
        lambda = Optim.minimizer(opts)

        #update internal branches
        for (i, node) in enumerate(nodes)
            if !isleaf(tree, node)
                tupdate = lambda[1] * t[i]
                getnodedata(tree, node).t = tupdate # tupdate is being put through as a vector
            end
        end


    end

    threepoint!(tree, trait, nodes)

    #information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q
    sigmahat = (nd.yy - 2 * betahat * nd.Q' + betahat * nd.xx * betahat') / n


    while any(i -> i < 0, sigmahat)
        leaves = getleaves(tree)
        for leaf in leaves
            ld = getnodedata(tree, leaf)
            ld.value = betahat - ld.value
        end

        threepoint!(tree, trait, nodes)
        sigmahat = nd.yy / n
    end

    
    negloglik = (1.0 / 2.0) * (n * log(2π) .+ nd.logV + n + n * log(det(sigmahat)))

    return betahat, sigmahat, negloglik, lambda #only return lambda if used
end

function estimaterates(tree::T, trait::Vector{String}; lambda = missing) where T <: AbstractTree
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
    NTraits = length(trait)

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

    return estimaterates!(tree, trait, lambda)
end

function tooptimise(lambda::Vector{Float64}, tree::T, nodes::Vector{N}, trait::Vector{String}, n::Int) where 
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B <: AbstractBranch{RT, NL},
    T <: AbstractTree{TT, RT, NL, N, B}}
    #lambda - value for Signal   
    #tree - tree
    #nodes - nodes of the tree in traversal order
    #t - vector of og branch lengths
    #trait - string of what trait is called on the tree
    #N - total number of nodes
    #n - number of leaves


    for node in nodes
    val = getnodedata(tree, node).value
    if hasinbound(tree, node)
    if isleaf(tree, node)
    #need parent to figure how much leaf branch needs to extend by - need special case if parent is the root node
    par = getparent(tree, node)

    #if par is root then len normal
    if isroot(tree, par)
        len = _getlength(tree, _getinbound(tree, node)) 
        td = traitdata(trait, val, len)
        setnodedata!(tree, node, td)

    else
        par_len = _getlength(tree, _getinbound(tree, par)) 

        #brach length is original leaf branch length, plus original parent branch length, minus new parent branch length
        len = _getlength(tree, _getinbound(tree, node)) + par_len - par_len * lambda[1] #get ansestor og branch length, add that then minus lambda*ansestor branch 
        td = traitdata(trait, val, len)
        setnodedata!(tree, node, td)
    end
    else
    len = _getlength(tree, _getinbound(tree, node)) * lambda[1]
    td = traitdata(trait, val, len)
    setnodedata!(tree, node, td)
    end
    else
    td = traitdata(trait, val)
    setnodedata!(tree, node, td)
    end
    end



    threepoint!(tree, trait, nodes)

    #information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q
    sigmahat = (nd.yy - 2 * betahat * nd.Q' + betahat * nd.xx * betahat') / n

    while any(i -> i < 0, sigmahat)
    leaves = getleaves(tree)
    for leaf in leaves
    ld = getnodedata(tree, leaf)
    ld.value = betahat - ld.value
    end

    threepoint!(tree, trait, nodes)
    sigmahat = nd.yy / n
    end

    negloglik = (1.0 / 2.0) * (n * log(2π) .+ nd.logV + n + n * log(det(sigmahat)))
    return negloglik
end

