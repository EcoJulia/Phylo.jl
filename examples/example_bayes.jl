using Turing
using Phylo
using Distributions
using Random
using BenchmarkTools
using DataFrames
using LinearAlgebra
using Plots 


#set seed
Random.seed!(678)

#generate a random tree
#number of tips on the tree
x=100

#generate tree
nu = Ultrametric{TraitTree{1}}(x);
tree = rand(nu);
distances(tree)

#get phylogentic variance matrix - needed for method using built in Julia functions
C = fill(1.0, (x,x)) - distances(tree)./2
C = abs.(Symmetric(C)) 

#generate traits, store in vector z
a = BrownianTrait(tree, "BMtrait", σ² = 1.0);
BMtraits = rand(a)

leafnames = getleafnames(tree);
z = Vector{Float64}();

for leaf in leafnames
    push!(z, BMtraits[leaf])
end

#using built in Julia functions
@model function phyloinf(z,C)
    beta ~ Uniform(-100, 100) 
    sigma ~ Uniform(0, 100)

    z ~ MvNormal(beta * ones(length(z)), sigma * C) 
end

c = sample(phyloinf(z, C), HMC(0.1, 5), 100000)

plot(c[:beta])
plot(c[:sigma])



#Use threepoint:
#pop trait data on tree
nodes = getnodes(tree, postorder)
rdat = DataFrame(species = leafnames, data = z)
for i in eachrow(rdat)
    setnodedata!(tree, i.species, Phylo.traitdata(["trait"], [i.data]));
end

#trait needs to be a vector of trait names, used for functions later
trait = ["trait"]


#add lengths to tree
for node in nodes
    val = getnodedata(tree, node).value
    if hasinbound(tree, node)
        len = Phylo.getlength(tree, Phylo.getinbound(tree, node)) 
        td = traitdata(trait, val, len)
        setnodedata!(tree, node, td)
    else
        td = traitdata(trait, val)
        setnodedata!(tree, node, td)
    end
end



@model function phyloinftree(tree, z) #z needs to be for leaves in postorder
    beta ~ Uniform(-100, 100) 
    sigma ~ Uniform(0, 100)
    
    z ~ Phylo.MyDist2(sigma, beta, tree) #tree.z ~ (impliment later)
end

c2 = sample(phyloinftree(tree, z), HMC(0.01, 5), 100000)
plot(c2[:beta])
plot(c2[:sigma])

#Signal
@model function phyloinftreelambda(tree, z) #z needs to be for leaves in postorder
    #calculate upper bound for signal
    intnodeheights = nodeheights(tree, noleaves=true)
    longnodeheight = maximum(intnodeheights)
    
    leafnodeheights = nodeheights(tree, onlyleaves=true)
    shortleafheight = minimum(leafnodeheights)
    
    upper = shortleafheight/longnodeheight

    lambda ~ Uniform(0, upper)
    beta ~ Uniform(-100, 100) 
    sigma ~ Uniform(0, 100)
    
    z ~ MyDist3(sigma, beta, lambda, tree) #tree.z ~ (impliment later)
end

c3 = sample(phyloinftreelambda(tree, z), HMC(0.01, 5), 10000) #add initial_params 
plot(c3[:beta])
plot(c3[:sigma])
plot(c3[:lambda])