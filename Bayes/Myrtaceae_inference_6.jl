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
using PDMats
using StatsPlots
using ReverseDiff

Random.seed!(1234)

# Bayesian inference on Myrtaceae data, 6 traits 

# Load the trees
include("read_tree_6.jl");


@model function βσ_multthreepoint(tree, z) 

    β ~ MvNormal(zeros(6), 1.0 * I)
    
    σ_scale ~ filldist(truncated(TDist(3); lower = 0), 6)  # positive half of t dist
    
    #Lcorr ~ LKJCholesky(6, 1.0)
    #σ = diagm(σ_scale) * Lcorr.L * Lcorr.L' * diagm(σ_scale)

    Lcorr ~ LKJ(6, 1.0)
    σ = Diagonal(σ_scale) * Lcorr * Diagonal(σ_scale)

    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end


n_samples = 100_000
n_warmup = 10_000


model6 = βσ_multthreepoint(tree6, z6);
#init = InitFromParams((β = zeros(6), σ_scale = ones(6), Lcorr = rand(LKJCholesky(6, 0.1))))
init = InitFromParams((β = zeros(6), σ_scale = ones(6), Lcorr = Matrix(I*1.0, 6, 6)))

chn6 = sample(model6, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; 
                initial_params = init,
                num_warmup = n_warmup)
# LKJCholesky 621.48 seconds
# LKJ 265.71 seconds <- before type changes, 621.61 seconds after

println(describe(chn6))
show(describe(chn6))
summary6 = DataFrame(chn6)
summary6desc = describe(summary6)

plot(chn6["β[1]"])
plot(chn6["σ_scale[1]"])
plot(chn6["Lcorr.L[6, 1]"])

est = estimaterates(tree6, trait)

sig_scale = summary6desc.mean[9:14]

v = summary6desc.mean[15:35]
corrl = zeros(6, 6)
corrl[tril!(trues(6, 6))] = v

sigma = diagm(sig_scale) * (corrl * corrl') * diagm(sig_scale)
est[2]

beta = summary6desc.mean[3:8]
est[1]


beta[1] * std(dat6.tmin) + mean(dat6.tmin)
beta[2] * std(dat6.tmax) + mean(dat6.tmax)
beta[3] * std(dat6.stl1) + mean(dat6.stl1)
beta[4] * std(dat6.swvl1) + mean(dat6.swvl1)
beta[5] * std(dat6.ssr) + mean(dat6.ssr)
beta[6] * std(dat6.tp) + mean(dat6.tp)

