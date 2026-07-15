# SPDX-License-Identifier: BSD-2-Clause

using DataFrames
using LinearAlgebra
using Optim
using Distributions
using ForwardDiff
using DynamicPPL

abstract type AbstractTraitData end

mutable struct TraitData{T <: Number, NTraits} <: AbstractTraitData # had to change types to Number for Bayes to work w/ lambda
    name::Vector{String}          # Name of traits
    value::Vector{Number}        # Trait values
    t::Number                      # branch length
    logV::Number
    p::Number                          # 1' V^(-1) 1
    yl::Vector{Number}                 # 1' V^(-1) y  (y is trait values)
    xl::Number                         # 1' V^(-1) x  (we set x as 1)
    Q::Vector{Number}                  # x' V^(-1) y
    xx::Number                         # x' V^(-1) x
    yy::Matrix{Number}                 # y' V^(-1) y
end

import Base.eltype
eltype(::TD) where {T <: Number, NT, TD <: TraitData{T, NT}} = T
eltype(::Type{TD}) where {T <: Number, NT, TD <: TraitData{T, NT}} = T

# change these from Float to <: Number too
function traitdata(::Type{T}, name::String, value::Float64,
                   t = 0.0) where {T <: Number}
    return traitdata(T, [name], [value], t)
end

function traitdata(::Type{T}, name::Vector{String}, value::AbstractArray,
                   t = 0.0) where {T}
    return TraitData{T, length(name)}(name, value, T(t), T(NaN), T(NaN),
                                      fill(T(NaN), length(name)), T(NaN),
                                      fill(T(NaN), length(name)), T(NaN),
                                      fill(NaN, length(name), length(name)), 
                                      )
end

function TraitData{T, NTraits}() where {T <: Number, NTraits}
    return traitdata(T, fill("", NTraits), fill(NaN, NTraits))
end

const TReB{T, NT, RT, BT, LenUnits} = RecursiveBranch{RT, String,
                                                      TraitData{T, NT}, Nothing,
                                                      BT, LenUnits}
const TReN{T, NT, RT, BT, LenUnits} = RecursiveNode{RT, String,
                                                    TraitData{T, NT}, Nothing,
                                                    BT, LenUnits}
const TReT{T, NT, RT, TD, BT, LenUnits} = RecursiveTree{RT, String,
                                                        TraitData{T, NT},
                                                        Nothing, BT, LenUnits,
                                                        TD}
const TReTD{T, NT, RT, BT, LenUnits} = TReT{T, NT, RT, Dict{String, Any}, BT,
                                            LenUnits}
const TraitTreeFloat64{NTraits} = TReTD{Float64, NTraits, OneRoot,
                                        PolytomousBranching, Float64}
const TraitTreeDual{NTraits} = TReTD{ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag,
                                                                      Float64},
                                                      Float64, 3}, NTraits,
                                     OneRoot, PolytomousBranching, Float64}
const TraitTreeNum{NTraits} = TReTD{Number, NTraits, OneRoot,
                                    PolytomousBranching, Float64}
const TraitTree{NTraits} = TReTD{Union{Float64,
                                       ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag,
                                                                        Float64},
                                                        Float64, 3}}, NTraits,
                                 OneRoot, PolytomousBranching, Float64}

function threepoint!(tree::T, trait::Vector{String},
                     nodes::Vector{N}) where
         {TT, RT, NL, N <: AbstractElt{RT, NL}, B <: AbstractElt{RT, NL},
          T <: AbstractTree{TT, RT, NL, N, B}}
    # prefrom algortihm in Ho & Ane 2014
    # function estimaterates gets inputs into right form
    # nodes - vector of nodes in the traversal order

    for node in nodes
        nd = getnodedata(tree, node)
        nodet = nd.t

        # need to see if node is a tip (leaf)
        if isleaf(tree, node)
            nodetrait = nd.value

            # update node data
            nd.logV = log(nodet)
            nd.p = inv(nodet)
            nd.yl = nodetrait
            nd.xl = 1.0
            nd.Q = nodetrait / nodet
            nd.xx = inv(nodet)
            nd.yy = (nodetrait * nodetrait') / nodet

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
            nd.yl = sum(ws .* [s.yl for s in childdata])

            c2 = nodet * pA^2 / calc
            nd.Q = sum(s.Q for s in childdata) - c2 * nd.xl * nd.yl
            nd.xx = sum(s.xx for s in childdata) - c2 * nd.xl * nd.xl
            nd.yy = sum(s.yy for s in childdata) - c2 * nd.yl * nd.yl'
        end
        # @assert getnodedata(tree, node) === nd
    end

    return tree
end

function estimaterates!(tree::T, trait::Vector{String},
                        lambda) where {T <: AbstractTree}

    # get information from tree in order to preform threepoint
    nodes = getnodes(tree, postorder)
    n = nleaves(tree)

    if lambda !== missing
        t = [getnodedata(tree, node).t for node in nodes]

        # info for optimiser
        lower = [0.0]
        start = [lambda]

        # upper is going to be largest lambda can be without going past the leaves
        # want to find longest internal branch
        intnodeheights = nodeheights(tree, noleaves = true)
        longnodeheight = maximum(intnodeheights)

        leafnodeheights = nodeheights(tree, onlyleaves = true)
        shortleafheight = minimum(leafnodeheights)

        upper = [shortleafheight / longnodeheight]

        # optimise to find lambda
        opts = optimize(x -> tooptimise(x, tree, nodes, trait, n),
                        lower, upper, start, Fminbox(LBFGS()),
                        Optim.Options(time_limit = 600))
        lambda = Optim.minimizer(opts)

        # update internal branches
        for (i, node) in enumerate(nodes)
            if !isleaf(tree, node)
                tupdate = lambda[1] * t[i]
                getnodedata(tree, node).t = tupdate
            end
        end
    end

    threepoint!(tree, trait, nodes)

    n = nleaves(tree)

    # information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q
    sigmahat = ((nd.yy .- 2 * nd.Q * betahat' .+ betahat * nd.xx * betahat') ./n)

    
    # NEED TO THINK ABOUT THIS
    #=
    while any(i -> i < 0, diag(sigmahat))
        leaves = getleaves(tree, postorder)
        for leaf in leaves
            ld = getnodedata(tree, leaf)
            ld.value = betahat - ld.value
        end

        threepoint!(tree, trait, nodes)
        sigmahat = nd.yy / n
    end
    =#
    k = length(trait)

    negloglik = (1.0 / 2.0) *
                (n * k * log(2π) + k * nd.logV + k * n + n * log(abs(det(sigmahat))))

    return betahat, sigmahat, negloglik, lambda # only return lambda if used
end

function estimaterates(tree::T, trait::Vector{String};
                       lambda = missing) where {T <: AbstractTree}
    # Returns evolution rate, starting value and negative log loglikelihood for traits on tip of tree
    # INPUTS
    # tree = tree with lengths, leaves all same length, trait data on leaves
    # trait = string with name of trait as found on leaves
    # OUTPUTS
    # sigmahat - evolution rate
    # betahat - estimated root trait value
    # negloglik - negative loglikelihood

    # need to add meaningful error when cant find traits on tree leaves

    nodes = getnodes(tree, anyorder)
    NTraits = length(trait)

    for node in nodes
        val = getnodedata(tree, node).value
        if hasinbound(tree, node)
            len = _getlength(tree, _getinbound(tree, node))
            td = traitdata(eltype(nodedatatype(T)), trait, val, len)
            setnodedata!(tree, node, td)
        else
            td = traitdata(eltype(nodedatatype(T)), trait, val)
            setnodedata!(tree, node, td)
        end
    end

    return estimaterates!(tree, trait, lambda)
end

function tooptimise(lambda::Vector{Float64}, tree::T, nodes::Vector{N},
                    trait::Vector{String},
                    n::Int) where
         {TT, RT, NL, N <: AbstractElt{RT, NL}, B <: AbstractElt{RT, NL},
          T <: AbstractTree{TT, RT, NL, N, B}}
    # lambda - value for Signal   
    # tree - tree
    # nodes - nodes of the tree in traversal order
    # t - vector of og branch lengths
    # trait - string of what trait is called on the tree
    # N - total number of nodes
    # n - number of leaves

    for node in nodes
        # val = getnodedata(tree, node).value
        if hasinbound(tree, node)
            len = _getlength(tree, _getinbound(tree, node)) * lambda[1]
            if isleaf(tree, node)
                full = getheight(tree, node)
                len = len + full * (1 - lambda[1])
                if len < 0
                    len = -len
                end
            end
            if len < 0
                len = -len
            end
            getnodedata(tree, node).t = len
        else
            getnodedata(tree, node).t = 0.0
        end
    end

    threepoint!(tree, trait, nodes)

    # information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q
    sigmahat = ((nd.yy .- 2 * nd.Q * betahat' .+ betahat * nd.xx * betahat') ./n)

    # with small numbers floating point errors may occur, if they do this should fix
    while any(i -> i < 0, sigmahat)
        leaves = getleaves(tree, postorder)
        for leaf in leaves
            ld = getnodedata(tree, leaf)
            ld.value = betahat - ld.value
        end

        threepoint!(tree, trait, nodes)
        sigmahat = nd.yy / n
    end
    

    negloglik = (1.0 / 2.0) *
                (n * log(2π) .+ nd.logV + n + n * log(det(sigmahat)))
    return negloglik
end

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

function threepointmultlambda!(tree::T, trait::Vector{String}, nodes::Vector{N},
                               C::Matrix{Float64},
                               lambda::Vector{Float64}) where
         {TT, RT, NL, N <: AbstractElt{RT, NL}, B <: AbstractElt{RT, NL},
          T <: AbstractTree{TT, RT, NL, N, B}}
    # prefrom algortihm in Ho & Ane 2014 for multiple lambda
    # function estimaterates gets inputs into right form

    for node in nodes
        nd = getnodedata(tree, node)
        nodet = nd.t

        C1 = diagm(sqrt.(lambda)) * C * diagm(sqrt.(lambda))

        C1inv = C1^(-1)

        k = length(trait)

        I = diagm(ones(k))

        if isleaf(tree, node)
            nodetrait = nd.value

            height = heighttoroot(tree, node)

            # update node data
            nd.logV = log(abs(det(C * height - (height - nodet) * C1)))
            nd.pmult = inv(height * C * C1inv - (height - nodet) * I)
            nd.ylmult = repeat(nodetrait, 1, k)
            nd.xlmult = ones(k)
            nd.Q = nd.ylmult' * nd.pmult * C1inv * ones(k)
            nd.xx = ones(k)' * nd.pmult * C1inv * ones(k)
            nd.yy = nd.ylmult' * nd.pmult * C1inv * nd.ylmult

        else
            # find direct desendents 
            children = getchildren(tree, node)

            # child data
            childdata = getnodedata.(tree, children)

            # calc pA and ws
            childp = [s.pmult for s in childdata]
            pA = sum(childp)

            nochildren = length(childp)

            ws = fill(fill(NaN, k, k), nochildren)

            pAinv = pA^(-1)
            for i in 1:nochildren
                ws[i] = pAinv * childp[i]
            end

            # update node data
            nd.logV = sum(s.logV for s in childdata) +
                      log(abs(det(I + nodet * pA))) # need to think about why the determinant is negative and if its okay
            nd.pmult = pA * (I + nodet * pA)^(-1)
            nd.xlmult = sum(ws .* [s.xlmult for s in childdata]) # dimensions may be wrong here
            nd.ylmult = sum(ws .* [s.ylmult for s in childdata])

            nd.Q = sum(s.Q for s in childdata) +
                   nd.ylmult' * nodet * pA^2 * (I + nodet * pA)^(-1) * C1inv *
                   nd.xlmult
            nd.xx = sum(s.xx for s in childdata) +
                    nd.xlmult' * (nodet * pA^2) * (I + nodet * pA)^(-1) *
                    C1inv * nd.xlmult
            nd.yy = sum(s.yy for s in childdata) +
                    nd.ylmult' * nodet * pA^2 * (I + nodet * pA)^(-1) * C1inv *
                    nd.ylmult
        end
        # @assert getnodedata(tree, node) === nd
    end

    return tree
end

function tooptimisemultlambda(lambda::Vector{Float64}, C::Matrix{Float64},
                              tree::T, nodes::Vector{N}, trait::Vector{String},
                              n::Int) where
         {TT, RT, NL, N <: AbstractElt{RT, NL}, B <: AbstractElt{RT, NL},
          T <: AbstractTree{TT, RT, NL, N, B}}
    # lambda - value for Signal 
    # C - evolutionary rates and covariances matrix  
    # tree - tree
    # nodes - nodes of the tree in traversal order
    # t - vector of og branch lengths
    # trait - string of what trait is called on the tree
    # N - total number of nodes
    # n - number of leaves

    threepointmultlambda!(tree, trait, nodes, C, lambda)
    # print after each call to double check runninng
    print('.')

    # information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q

    k = length(trait)

    negloglik = (1.0 / 2.0) * (n * k * log(2π) + nd.logV + n +
                 log(abs(det((nd.yy - 2 * betahat * nd.Q' +
                              betahat * nd.xx * betahat')))))
    return negloglik
end

function estimaterates!(tree::T, trait::Vector{String},
                        lambda::Vector{Float64}) where {T <: AbstractTree}

    # gather information from tree in order to preform threepoint
    nodes = getnodes(tree, postorder)
    n = nleaves(tree)

    # Get C by running algorithm w/o lambda
    betahat, C, negloglik = estimaterates(tree, trait)
    print(C)

    t = [getnodedata(tree, node).t for node in nodes]

    k = length(lambda)

    # info for optimiser
    lower = fill(floatmin(), k)

    # upper is going to be largest lambda can be without going past the leaves
    # want to find longest internal branch
    intnodeheights = nodeheights(tree, noleaves = true)
    longnodeheight = maximum(intnodeheights)

    leafnodeheights = nodeheights(tree, onlyleaves = true)
    shortleafheight = minimum(leafnodeheights)

    upper = fill(shortleafheight / longnodeheight, k)

    # lowest value of C, diagonals cant be lower than zero, off diagonals dont have restriction 
    lowerC = fill(-Inf, k, k)

    for i in 1:k
        lowerC[i, i] = 0
    end

    # upper value of C, no restirction
    upperC = fill(Inf, k, k)

    for i in 1:3
        # optimise to find lambda
        optslambda = optimize(x -> tooptimisemultlambda(x, C, tree, nodes,
                                                        trait, n),
                              lower, upper, lambda, Fminbox(LBFGS()),
                              Optim.Options(time_limit = 30))
        # get lambda value
        lambda = Optim.minimizer(optslambda)
        # print to ensure its running and see how it's going
        print(lambda)

        # optimise to find C
        optsC = optimize(x -> tooptimisemultlambda(lambda, x, tree, nodes,
                                                   trait, n),
                         lowerC, upperC, C, Fminbox(LBFGS()),
                         Optim.Options(time_limit = 30))
        # get C value
        C = Optim.minimizer(optsC)
    end

    # information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q

    k = length(trait)

    negloglik = (1.0 / 2.0) * (n * k * log(2π) + nd.logV + n +
                 log(abs(det((nd.yy - 2 * betahat * nd.Q' +
                              betahat * nd.xx * betahat')))))

    return betahat, C, negloglik, lambda
end



