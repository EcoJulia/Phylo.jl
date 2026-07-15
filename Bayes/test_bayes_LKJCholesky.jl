using Phylo
using DataFrames
using ForwardDiff
using DynamicPPL
using Turing
using LinearAlgebra
using Plots
using PDMats
using Random

# Bayesian analysis test using generated data

Random.seed!(1234)

# set number of tips of the tree
n_tips = 100

# generate a tree
nu = Ultrametric{TraitTree{3}}(n_tips)
tree = rand(nu);

σ = [1.0 0.0 0.0;
    0.0 1.0 0.0;
    0.0 0.0 1.0]

a = BrownianTraitMult(tree, ["trait1", "trait2", "trait3"], [2.0, -2.0, 1.0]; σ = σ);
traits = rand(a);


leafnames = getleafnames(tree);
# add data to the tree
for i in leafnames
    setnodedata!(tree, i, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait1", "trait2", "trait3"], traits[i]))
end

# add lengths to tree
trait = ["trait1", "trait2", "trait3"]
nodes = getnodes(tree, postorder)

for node in nodes
    val = getnodedata(tree, node).value
    if hasinbound(tree, node)
        len = Phylo.getlength(tree, Phylo.getinbound(tree, node))
        td = traitdata(eltype(nodedatatype(typeof(tree))), trait, val, len)
        setnodedata!(tree, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(tree))), trait, val)
        setnodedata!(tree, node, td)
    end
end

# get z in the form we want
zvec = [getnodedata(tree, leaf).value for leaf in leafnames]
z = reduce(vcat, zvec)


@model function βσ_multthreepoint(tree, z) 

    β ~ MvNormal(zeros(3), 1.0 * I)
    
    σ_scale ~ filldist(truncated(TDist(3); lower = 0), 3)  # positive half of t dist
    Lcorr ~ LKJCholesky(3, 1.0)
    Σ_L = diagm(σ_scale) * Lcorr.L
    σ = Σ_L * Σ_L'

    z ~ Phylo.MyDist4(σ, β, tree) 

    return nothing
end


n_samples = 100_000
n_warmup = 10_000


model = βσ_multthreepoint(tree, z);

init = InitFromParams((β = zeros(3), σ_scale = ones(3), Lcorr = rand(LKJCholesky(3, 0.1))))

chn = sample(model, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; 
                initial_params=init,
                num_warmup = n_warmup)
# 58.13 secs
println(describe(chn))

summary6_1 = describe(chn)
summary_stats6_1 = summary6_1[1]
summarydf6_1 = DataFrame(summary6_1)
println(summarydf6_1)
println(DataFrame(summary6_1[2]))

plot(chn["β[1]"])
plot(chn["σ_scale[1]"])
plot(chn["Lcorr.L[3, 1]"])
plot(chn["Lcorr.L[2, 2]"])
plot(chn["Lcorr.L[3, 2]"])
plot(chn["Lcorr.L[3, 3]"])

beta = mean(group(chn, :β))[:,2]

sigscale = mean(group(chn, :σ_scale))[:,2]

v = mean(group(chn, "Lcorr.L"))[:,2]
corrl = zeros(3, 3)
corrl[tril!(trues(3, 3))] = v

sigma = diagm(sigscale) * (corrl * corrl') * diagm(sigscale)