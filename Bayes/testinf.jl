using Turing
using Phylo
using Distributions
using Random
using BenchmarkTools
using DataFrames
using LinearAlgebra
using Plots 

#set seed
Random.seed!(18)

x=10
#generate tree
nu = Ultrametric{TraitTree{1}}(x);
tree = rand(nu)

#get phylogentic variance matrix
C = fill(2.0, (x,x)) - distances(tree)
C = Symmetric(C)


a = BrownianTrait(tree, "BMtrait")
BMtraits = rand(a)

leafnames = getleafnames(tree);
z = Vector{Float64}();

for leaf in leafnames
    push!(z, BMtraits[leaf])
end

@model function phyloinf(z,C)
    x ~ Uniform(-100, 100) 
    sigma ~ Uniform(0, 100)

    z ~ MvNormal(x * ones(length(z)), sigma * C)  #likelihood calculated using logpdf
    return x, sigma
end

c = sample(phyloinf(z, C), HMC(0.1, 5), 1000)


#use my code to calculate loglikelihood
#we know sigma and n so just need to calculate logV <- will only need to do once in non signal inference
nodes = getnodes(tree, postorder)

#Crate dataframe
rdat = DataFrame(species = leafnames, data = z)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree, i.species, Phylo.traitdata(["trait"], [i.data]));
end

#get nodes
nodes = getnodes(tree, postorder)

#add lengths to tree
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

n = nleaves(tree)

@model function phyloinftree(tree, trait, nodes, n)
    x ~ Uniform(-100, 100) 
    sigma ~ Uniform(0, 100)

    threepoint!(tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(tree, nN)

    negloglik = (1.0 / 2.0) * (n * log(2π) .+ nd.logV + n + n * log(det(sigma))) #need beta/x in this
    return x, sigma
end
