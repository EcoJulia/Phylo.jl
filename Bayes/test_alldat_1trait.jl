using Phylo
using DataFrames
using ForwardDiff
using DynamicPPL
using Turing
using LinearAlgebra
using Plots
using PDMats
using Random

# comparisons on methods to include all data


Random.seed!(1234)

# set number of tips of the tree
n_tips = 100

# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
tree = rand(nu);

σ = 1.0

a = BrownianTrait(tree, "trait", 10.0; σ = σ);
traits = rand(a);


leafnames = getleafnames(tree);

# add lengths to tree
trait = ["trait"]
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

leaves = getleaves(tree, postorder) #needs to be postorder so in the same form as z
# create dataframe for values

# number of values at each tip
yn = 100

df = DataFrame(leaf = repeat(leaves, inner = yn), 
               value= zeros(length(leaves)*yn))

sigma_y = 0.1

for i in eachrow(df)
    leafname = getnodename(tree,i[1])
    y = rand(Normal(traits[leafname], sigma_y))
    i[2] = y
end

df

transform!(groupby(df, :leaf), groupindices => :leafnumber)

number = df[!, :leafnumber]
values = df[!, :value]

@model function alldata(tree, value, number, n_tips) 

    β ~ Normal(0.0, 10.0)
    
    σ ~ truncated(TDist(3); lower = 0) # positive half of t dist

    z ~ Phylo.MyDist2(σ, β, tree) 

    σ_y ~ filldist(truncated(TDist(3); lower = 0), n_tips) # positive half of t dist
    zall = z[number]
    σ_yall = σ_y[number]
    value ~ MvNormal(zall, I*σ_yall)
    

    return nothing
end

#Turing.@addlogprob! sum(logpdf.(Normal.(μ, σ), y))

n_samples = 100_000
n_warmup = 10_000


model = alldata(tree, values, number, n_tips);
init = InitFromParams((β = 0.0, σ = 1.0, z = fill(0.0,n_tips), σ_y = fill(1.0,n_tips)))
# ; adtype=AutoReverseDiff(; compile=true)
chn = sample(model, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; initial_params=init, num_warmup = n_warmup)
# 387.86 seconds


println(describe(chn))
show(describe(chn))
summary = DataFrame(chn)
summarydesc = describe(summary)

plot(chn["β"])
plot(chn["σ"])
plot(chn["z[1]"])
plot(chn["z[2]"])
plot(chn["σ_y[1]"])
plot(chn["σ_y[2]"])

estz = summarydesc.mean[5:14]
truez = [traits[leaf] for leaf in getleafnames(tree, postorder)]



# using sufficient statistics

# number of values for each tip
n_values = [count(==(element),df.leafnumber) for element in unique(df.leafnumber)]

# sum of values for each tip
sum_values = zeros(n_tips)

for i in 1:n_tips
    sum_values[i] = sum(df[df.leafnumber .== i, :value])
end

sum_values

# sum of squares for each tip
sum_sq_values = zeros(n_tips)

for i in 1:n_tips
    sum_sq_values[i] = sum(df[df.leafnumber .== i, :value].^2)
end

sum_sq_values

#@addlogprob! 

z = rand(n_tips)

σ_y = rand(filldist(truncated(TDist(3); lower = 0), n_tips))

sum( @. (-(n_values/2) * log(2π) - n_values * log(σ_y) - (sum_sq_values - 2 * z * sum_values + n_values * z^2) / (2 * σ_y^2)))

@model function alldata2(tree, n_values, sum_sq_values, sum_values, n_tips) 

    β ~ Normal(0.0, 10.0)
    
    σ ~ truncated(TDist(3); lower = 0) # positive half of t dist

    z ~ Phylo.MyDist2(σ, β, tree) 

    σ_y ~ filldist(truncated(TDist(3); lower = 0), n_tips) # positive half of t dist
    
    # value ~ MvNormal(zall, I*σ_yall)
    @addlogprob! sum( @. (-(n_values/2) * log(2π) - n_values * log(σ_y) 
                      - (sum_sq_values - 2 * z * sum_values + n_values * z^2) / (2 * σ_y^2)))
    

    return nothing
end

#Turing.@addlogprob! sum(logpdf.(Normal.(μ, σ), y))

n_samples = 100_000
n_warmup = 10_000


model2 = alldata2(tree, n_values, sum_sq_values, sum_values, n_tips); 
init = InitFromParams((β = 0.0, σ = 1.0, z = fill(0.0,n_tips), σ_y = fill(1.0,n_tips)))
# ; adtype=AutoReverseDiff(; compile=true)
chn2 = sample(model2, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; initial_params=init, num_warmup = n_warmup)
# 10 tips 37.34 seconds
# 100 tips 332.17 seconds

println(describe(chn2))
show(describe(chn2))
summary2 = DataFrame(chn2)
summarydesc2 = describe(summary2)

# try with mean and ssd (sum squared deviation)

# number of values for each tip
n_values = [count(==(element),df.leafnumber) for element in unique(df.leafnumber)]

# means
mean_values = zeros(n_tips)

for i in 1:n_tips
    mean_values[i] = mean(df[df.leafnumber .== i, :value])
end

mean_values

# ssd for each tip
ssd_values = zeros(n_tips)

for i in 1:n_tips
    ssd_values[i] = sum((df[df.leafnumber .== i, :value] .- mean_values[i]).^2)
end

ssd_values

# test log likelihood
#z = rand(10)

#σ_y = rand(filldist(truncated(TDist(3); lower = 0), 10))

sum( @. (-(n_values/2) * log(2π) - n_values * log(σ_y) - (ssd_values + n_values * (mean_values - z)^2) / (2 * σ_y^2)))

@model function alldata3(tree, n_values, mean_values, ssd_values, n_tips) 

    β ~ Normal(0.0, 10.0)
    
    σ ~ truncated(TDist(3); lower = 0) # positive half of t dist

    z ~ Phylo.MyDist2(σ, β, tree) 

    σ_y ~ filldist(truncated(TDist(3); lower = 0), n_tips) # positive half of t dist
    
    # value ~ MvNormal(zall, I*σ_yall)
    @addlogprob! sum( @. (-(n_values/2) * log(2π) - n_values * log(σ_y) 
                      - (ssd_values + n_values * (mean_values - z)^2) / (2 * σ_y^2)))
    

    return nothing
end

n_samples = 100_000
n_warmup = 10_000


model3 = alldata3(tree, n_values, mean_values, ssd_values, n_tips); 
init = InitFromParams((β = 0.0, σ = 1.0, z = fill(0.0,n_tips), σ_y = fill(1.0,n_tips)))
# ; adtype=AutoReverseDiff(; compile=true)
chn3 = sample(model3, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; initial_params=init, num_warmup = n_warmup)
# 10 tips 62.61 seconds
# 100 tips 309.81 seconds

println(describe(chn3))
show(describe(chn3))
summary2 = DataFrame(chn3)
summarydesc2 = describe(summary3)