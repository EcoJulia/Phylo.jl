using DataFrames
using LinearAlgebra
using Optim
using StaticArrays


abstract type AbstractTraitData end

mutable struct TraitData{NTraits} <: AbstractTraitData 
    name::Vector{String}          #Name of traits
    value::Vector{Float64}        #Trait values
    t::Float64                    #branch length
    logV::Float64
    p::Float64                    # 1' V^(-1) 1
    yl::Vector{Float64}           # 1' V^(-1) y  (y is trait values)
    xl::Float64                   # 1' V^(-1) x  (we set x as 1)
    Q::Vector{Float64}            # x' V^(-1) y
    xx::Float64                   # x' V^(-1) x
    yy::Matrix{Float64}           # y' V^(-1) y
    v::Vector{Float64}            #used for PIC
    xlmult::Vector{Float64}       # 1' W^(-1) x  
    pmult::Matrix{Float64}        # 1' W^(-1) C1      
    xxmult::Matrix{Float64}       # x' W^(-1) x 
    yymult::Float64               # y' W^(-1) y
    Qmult::Float64                # x' W^(-1) y
    ylmult::Matrix{Float64}       # 1' W^(-1) y
end

traitdata(name::Vector{String}, value::Vector{Float64}, t::Float64 = 0.0) =
    TraitData{length(name)}(name, value, t, NaN, NaN, fill(NaN, length(name)), NaN, 
                            fill(NaN, length(name)), NaN, fill(NaN, length(name), length(name)), fill(NaN, length(name)), fill(NaN, length(name)),
                            fill(NaN, length(name), length(name)), fill(NaN, length(name), length(name)), NaN, NaN, fill(NaN, length(name), length(name)))

TraitData{NTraits}() where NTraits = traitdata(fill("", NTraits), fill(NaN, NTraits))

const TLB{RT, LenUnits} = LinkBranch{RT, String, Nothing, LenUnits}
const TLN{NT, RT, LenUnits} = LinkNode{RT, String, TraitData{NT}, TLB{RT, LenUnits}}
const TLT{NT, RT, TD, LenUnits} = LinkTree{RT, String, TLN{NT, RT, LenUnits}, TLB{RT, LenUnits}, TD}
const TLTD{NT, RT, LenUnits} = TLT{NT, RT, Dict{String, Any}, LenUnits}
const TraitTree{NTraits} = TLTD{NTraits, OneRoot, Float64}



function threepoint!(tree::T, trait::Vector{String}, nodes::Vector{N}) where 
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

    return tree
end

function estimaterates!(tree::T, trait::Vector{String}, lambda) where T <: AbstractTree

    #get information from tree in order to preform threepoint
    nodes = getnodes(tree, postorder)
    

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
                getnodedata(tree, node).t = tupdate 
            end
        end


    end

    threepoint!(tree, trait, nodes)

    n = nleaves(tree)

    #information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q
    sigmahat = ((nd.yy - 2 * betahat * nd.Q' + betahat * nd.xx * betahat') / n)

 #NEED TO THINK ABOUT THIS
    while any(i -> i < 0, diag(sigmahat))
        leaves = getleaves(tree)
        for leaf in leaves
            ld = getnodedata(tree, leaf)
            ld.value = betahat - ld.value
        end

        threepoint!(tree, trait, nodes)
        sigmahat = nd.yy / n
    end

    k = length(trait)

    negloglik = (1.0 / 2.0) * (n * k * log(2π) + nd.logV + n + n * log(abs(det(sigmahat)))) 


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
        #val = getnodedata(tree, node).value
        if hasinbound(tree, node)
            len = _getlength(tree, _getinbound(tree, node)) * lambda[1]
            if isleaf(tree, node)
                full = getheight(tree, node)
                len = len + full * (1 - lambda[1])
            end
            getnodedata(tree, node).t = len
        else
            getnodedata(tree, node).t = 0
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

####################################################################

function threepointmultlambda!(tree::T, trait::Vector{String}, nodes::Vector{N}, C::Matrix{Float64}, lambda::Vector{Float64}) where 
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}}
    #prefrom algortihm in Ho & Ane 2014 for multiple lambda
    #function estimaterates gets inputs into right form

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
            nd.ylmult = repeat(nodetrait, 1, k)'
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
            nd.logV = sum(s.logV for s in childdata) + log(abs(det(I + nodet * pA))) #need to think about why the determinant is negative and if its okay
            nd.pmult = pA * (I + nodet * pA)^(-1)
            nd.xlmult = sum(ws .* [s.xlmult for s in childdata]) #dimensions may be wrong here
            nd.ylmult = sum(ws .* [s.ylmult for s in childdata])


            nd.Q = sum(s.Q for s in childdata) + nd.ylmult' * nodet * pA^2 * (I + nodet * pA)^(-1) * C1inv * nd.xlmult
            nd.xx = sum(s.xx for s in childdata) + nd.xlmult' * (nodet * pA^2) * (I + nodet * pA)^(-1) * C1inv * nd.xlmult
            nd.yy = sum(s.yy for s in childdata) + nd.ylmult' * nodet * pA^2 * (I + nodet * pA)^(-1) * C1inv * nd.ylmult


        end
        # @assert getnodedata(tree, node) === nd
    end

    return tree
end

function tooptimisemultlambda(lambda::Vector{Float64}, C::Matrix{Float64}, tree::T, nodes::Vector{N}, trait::Vector{String}, n::Int) where 
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B <: AbstractBranch{RT, NL},
    T <: AbstractTree{TT, RT, NL, N, B}}
    #lambda - value for Signal 
    #C - evolutionary rates and covariances matrix  
    #tree - tree
    #nodes - nodes of the tree in traversal order
    #t - vector of og branch lengths
    #trait - string of what trait is called on the tree
    #N - total number of nodes
    #n - number of leaves


    threepointmultlambda!(tree, trait, nodes, C, lambda)
    #print after each call to double check runninng
    print('.')


    #information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q

    k = length(trait)

    negloglik = (1.0 / 2.0) * (n * k * log(2π) + nd.logV + n + log(abs(det((nd.yy - 2 * betahat * nd.Q' + betahat * nd.xx * betahat')))))
    return negloglik
end

function estimaterates!(tree::T, trait::Vector{String}, lambda::Vector{Float64}) where T <: AbstractTree

    #gather information from tree in order to preform threepoint
    nodes = getnodes(tree, postorder)
    n = nleaves(tree)

    #Get C by running algorithm w/o lambda
    betahat, C, negloglik = estimaterates(tree, trait) 
    print(C)
    
    t = [getnodedata(tree, node).t for node in nodes]

    k = length(lambda)

    #info for optimiser
    lower = fill(floatmin(), k)

    #upper is going to be largest lambda can be without going past the leaves
    #want to find longest internal branch
    intnodeheights = nodeheights(tree, noleaves=true)
    longnodeheight = maximum(intnodeheights)

    leafnodeheights = nodeheights(tree, onlyleaves=true)
    shortleafheight = minimum(leafnodeheights)

    upper = fill(shortleafheight/longnodeheight, k)
    
    #lowest value of C, diagonals cant be lower than zero, off diagonals dont have restriction 
    lowerC = fill(-Inf, k, k)

    for i in 1:k
        lowerC[i,i] = 0
    end

    #upper value of C, no restirction
    upperC = fill(Inf, k, k)


    for i in 1:3
        #optimise to find lambda
        optslambda = optimize(x -> tooptimisemultlambda(x, C, tree, nodes, trait, n),
                   lower, upper, lambda, Fminbox(LBFGS()), Optim.Options(time_limit = 30))
        #get lambda value
        lambda = Optim.minimizer(optslambda)
        #print to ensure its running and see how it's going
        print(lambda)

        #optimise to find C
        optsC = optimize(x -> tooptimisemultlambda(lambda, x, tree, nodes, trait, n), 
                    lowerC, upperC, C, Fminbox(LBFGS()), Optim.Options(time_limit = 30))
        #get C value
        C = Optim.minimizer(optsC)
        
    end


    #information from last node
    nN = last(nodes)
    nd = getnodedata(tree, nN)

    betahat = inv(nd.xx) * nd.Q

    k = length(trait)

    negloglik = (1.0 / 2.0) * (n * k * log(2π) + nd.logV + n + log(abs(det((nd.yy - 2 * betahat * nd.Q' + betahat * nd.xx * betahat')))))

    return betahat, C, negloglik, lambda
end


#PIC calculations, ignore for now
##############################
#=
function pic!(tree::T, trait::Vector{String}, nodes::Vector{N}, lambda::Vector{Float64}) where 
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}}
    #prefrom algortihm in Ho & Ane 2014
    #function estimaterates gets inputs into right form
    #nodes - vector of nodes in the traversal order

    #number of internal nodes
    internal = nleaves(tree)

    #number or traits
    k = length(trait)

    #get length of whole tree
    full = getheight(tree, node)

    #create empty vectors, u and V to store things in for later
    U = fill(fill(NaN, k), internal)
    V = fill(fill(NaN, k), internal)

    #use i to count which elements of u and V to be updated
    i=1

    for node in nodes
        nd = getnodedata(tree, node)
        nodet = nd.t

        #need to see if node is a tip (leaf)
        if isleaf(tree, node)

            #for Freckleton likelihood calculations
            nd.v = nodet * ones(k) .+ full * ones(k) .* (ones(k) .- lambda)

        else
            # need to find direct desendents 
            children = getchildren(tree, node)

            # child data
            childdata = getnodedata.(tree, children)

            #for Freckleton likelihood calculations
            #assuming two children
            child1value = childdata[1].value
            child2value = childdata[2].value
            child1v = childdata[1].v
            child2v = childdata[2].v

            U[i] = abs.(child1value - child2value) 
            V[i] = child1v .+ child2v

            nd.value = ((child1value ./ child1v) .+ (child2value ./ child2v)) ./ ((1 / child1v) .+ (1 / child2v))
            nd.v = nodet * ones(k) .* lambda .+ (child1v .* child2v) ./ (child1v .+ child2v)

            i=i+1
        end
        # @assert getnodedata(tree, node) === nd
    end
end

    function estimateratespic!(tree::T, trait::Vector{String}, lambda) where T <: AbstractTree

        #get information from tree in order to preform threepoint
        nodes = getnodes(tree, postorder)
        n = nleaves(tree)
    
        a = threepointflambda!(tree, trait, nodes)
    
        tree = a[1]
        U = a[2]
        V = a[3]
    
        #information from last node
        nN = last(nodes)
        nd = getnodedata(tree, nN)
    
        #Freckleton calculations
        k = length(trait)
    
        UU = fill(fill(NaN, k, k), n)
    
        for i in 1:n
            UU[i] = U[i] * U[i]' / V[i]
        end
    
        sigma2 = 1/(n) * sum(UU) #I have wrong off diagonals  - try and use R code?
    
    
        UCU = fill(NaN, n)
    
       for i in 1:n
            UCU[i] = ((V[i] .* U[i])' * sigma2^(-1) * (V[i] .* U[i]))  
        end
    
        nll2 = (1.0 / 2.0) * (n * k * log(2π) + n * log(det(sigma2)) + sum(log.(V)) + sum(UCU)) 
    
        beta = nd.value
    
        return lambda, beta, sigma2, nll2 #only return lambda if used
    end
    
    function estimateratespic(tree::T, trait::Vector{String}; lambda = missing) where T <: AbstractTree
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

        if lambda == missing
            lambda = ones(NTraits)
        end
    
        return estimateratesflambda!(tree, trait, lambda)
    end

=#
    
    