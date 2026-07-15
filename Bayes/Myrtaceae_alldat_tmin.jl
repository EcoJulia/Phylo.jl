# SPDX-License-Identifier: BSD-2-Clause

using Turing
using Phylo
using Distributions
using Random
using DataFrames
using LinearAlgebra
using Plots
using CSV
using ForwardDiff
using DynamicPPL
using StatsPlots
using ReverseDiff

# Bayesian analysis on Myrtaceae tree using all data avalible for tmin

# Load the trees
include("read_tree_alldat_tmin.jl");

n_leaves = nleaves(tree)

# number of values for each tip
n_values = [count(==(element),df.leafnumber) for element in unique(df.leafnumber)]

# means
mean_values = zeros(n_leaves)

for i in 1:n_leaves
    mean_values[i] = mean(df[df.leafnumber .== i, :tmin])
end

mean_values

# ssd for each tip
ssd_values = zeros(n_leaves)

for i in 1:n_leaves
    ssd_values[i] = sum((df[df.leafnumber .== i, :tmin] .- mean_values[i]).^2)
end

ssd_values


@model function alldata1(tree, n_values, mean_values, ssd_values) 

    β ~ Normal(285.0, 10.0)
    
    σ ~ truncated(TDist(10); lower = 0) # positive half of t dist

    z ~ Phylo.MyDist2(σ, β, tree) 

    σ_y ~ truncated(TDist(3); lower = 0)
    

    # values ~ MvNormal(zall, I*σ_y)
    @addlogprob! sum( @. (-(n_values/2) * log(2π) - n_values * log(σ_y) 
                          - (ssd_values + n_values * (mean_values - z)^2) / (2 * σ_y^2)))

    return nothing
end

n_samples = 10_000
n_warmup = 1_000


model1 = alldata1(tree, n_values, mean_values, ssd_values);
#model1 = alldata1(tree, df);
init1 = InitFromParams((β = 285.0, σ = 5.0, z = fill(285.0,n_leaves)))
# ; adtype=AutoReverseDiff(; compile=true)
chn1 = sample(model1, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; num_warmup = n_warmup)
# 100_000 samples
# reduced time 1682.02 seconds ~ 30 mins
# alldat (610,726 entries) time 170 hours
# using means and ssd 24 hours
# 10_000 samples 5427.04 seconds ~ 1.5 hours
println(describe(chn1))
summary1 = DataFrame(chn1)
summarydesc1 = describe(summary1)

summarydesc1[summarydesc1.variable .== :σ_y, :]

plot(chn1["β"])
plot(chn1["σ"])
plot(chn1["z[1]"])
plot(chn1["z[2]"])
plot(chn1["σ_y"])

sigma_y = summarydesc1.mean[544]

@model function alldata2(tree, n_leaves, n_values, mean_values, ssd_values) 

    β ~ Normal(0.0, 10.0)
    
    σ ~ truncated(TDist(3); lower = 0) # positive half of t dist

    z ~ Phylo.MyDist2(σ, β, tree) 

    σ_y ~ filldist(truncated(TDist(3); lower = 0), n_leaves) # positive half of t dist

    #df[!, :tmin] ~ MvNormal(zall, I*σ_yall)
    @addlogprob! sum( @. (-(n_values/2) * log(2π) - n_values * log(σ_y) 
                          - (ssd_values + n_values * (mean_values - z)^2) / (2 * σ_y^2)))

    return nothing
end


n_samples = 20_000
n_warmup = 5_000


model2 = alldata2(tree, n_leaves, n_values, mean_values, ssd_values); 
#model2 = alldata1(tree, df, numleaves);
init2 = InitFromParams((β = 285.0, σ = 5.0, z = fill(285.0,n_leaves)))
# ; adtype=AutoReverseDiff(; compile=true)
chn2 = sample(model2, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; initial_params=init2, num_warmup = n_warmup)
# 7080.61 seconds ~ 2 hours
println(describe(chn2))
show(describe(chn2))
summary2 = DataFrame(chn2)
summarydesc2 = describe(summary2)

plot(chn2["β"])
plot(chn2["σ"])
plot(chn2["z[1]"])
plot(chn2["z[2]"])
plot(chn2["σ_y[1]"])
plot(chn2["σ_y[2]"])

