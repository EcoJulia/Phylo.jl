using Phylo
using DataFrames
using Random
using CSV
using Optim #need to remember to remove

estimaterates(multtree, ["trait1", "trait2"], lambda = [0.1, 0.1])

Random.seed!(18)

x=10
#generate tree
nu = Ultrametric{TraitTree{2}}(x);
tree = rand(nu)

a = BrownianTrait(tree, "BMtrait")
b = BrownianTrait(tree, "BMtrait")
BMtraits1 = rand(a)
BMtraits2 = rand(b)

leafnames = getleafnames(tree);
z1 = Vector{Float64}();
z2 = Vector{Float64}();

for leaf in leafnames
    push!(z1, BMtraits1[leaf])
    push!(z2, BMtraits2[leaf])
end

nodes = getnodes(tree, postorder)

#Crate dataframe
rdat = DataFrame(species = leafnames, data1 = z1, data2 = z2)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree, i.species, Phylo.traitdata(["trait1", "trait2"], [i.data1, i.data2]));
end

betahat, C, negloglik = estimaterates(tree, ["trait1", "trait2"]) 

nodes = getnodes(tree)

trait = ["trait1", "trait2"]

n = nleaves(tree)

lambda = [0.8, 0.8]

lower = zeros(length(lambda))
    
#upper is going to be largest lambda can be without going past the leaves
#want to find longest internal branch
intnodeheights = nodeheights(tree, noleaves=true)
longnodeheight = maximum(intnodeheights)

leafnodeheights = nodeheights(tree, onlyleaves=true)
shortleafheight = minimum(leafnodeheights)

upper = fill(shortleafheight/longnodeheight, length(lambda))

optslambda = optimize(x -> tooptimisemultlambda(x, C, tree, nodes, trait, n),
                       lower, upper, lambda, Fminbox(LBFGS()), Optim.Options(iterations = 300))

lambda = Optim.minimizer(optslambda)





Random.seed!(18)

x=10
#generate tree
nu = Ultrametric{TraitTree{2}}(x);
tree = rand(nu)

a = BrownianTrait(tree, "BMtrait")
b = BrownianTrait(tree, "BMtrait")
BMtraits1 = rand(a)
BMtraits2 = rand(b)

leafnames = getleafnames(tree);
z1 = Vector{Float64}();
z2 = Vector{Float64}();

for leaf in leafnames
    push!(z1, BMtraits1[leaf])
    push!(z2, BMtraits2[leaf])
end

nodes = getnodes(tree, postorder)

#Crate dataframe
rdat = DataFrame(species = leafnames, data1 = z1, data2 = z2)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree, i.species, Phylo.traitdata(["trait1", "trait2"], [i.data1, i.data2]));
end

estimaterates(tree, ["trait1", "trait2"], lambda = [0.8, 0.8])