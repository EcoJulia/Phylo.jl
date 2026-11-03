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
include("read_tree_6_nontransformed.jl");


@model function βσ_multthreepoint(tree, z) 

    sigmabeta = diagm([8.0, 8.0, 8.0, 0.1, 2.0e6, 0.01])
    meanbeta = [280.0, 290.0, 290.0, 0.1, 1.6e7, 0.01]
    β ~ MvNormal(meanbeta, sigmabeta)
    
    σ_scale ~ filldist(truncated(TDist(3); lower = 0), 6)  # positive half of t dist
    Lcorr ~ LKJCholesky(6, 1.0)
    σ = diagm(σ_scale) * Lcorr.L * Lcorr.L' * diagm(σ_scale)

    #Lcorr ~ LKJ(6, 1.0)
    #σ = Diagonal(σ_scale) * Lcorr * Diagonal(σ_scale)

    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end


n_samples = 100_000
n_warmup = 10_000

beta = [mean(dat6.tmin), mean(dat6.tmax), mean(dat6.stl1), mean(dat6.swvl1), mean(dat6.ssr), mean(dat6.tp)]
model6 = βσ_multthreepoint(tree6nt, z6nt);
init = InitFromParams((β = beta, σ_scale = ones(6), Lcorr = rand(LKJCholesky(6, 0.1))))
#init = InitFromParams((β = zeros(6), σ_scale = ones(6), Lcorr = Matrix(I*1.0, 6, 6)))

chn6nt = sample(model6, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; 
                initial_params = init,
                num_warmup = n_warmup)
println(describe(chn6nt))
summary6nt = DataFrame(chn6nt)
summary6ntdesc = describe(summary6nt)

plot(chn6nt["β[1]"])
plot(chn6nt["σ_scale[1]"])
plot(chn6nt["Lcorr.L[6, 1]"])

estnt = estimaterates(tree6nt, trait)

sig_scalent = summary6ntdesc.mean[9:14]

vnt = summary6ntdesc.mean[15:35]
corrlnt = zeros(6, 6)
corrlnt[tril!(trues(6, 6))] = vnt

sigmant = diagm(sig_scalent) * (corrlnt * corrlnt') * diagm(sig_scalent)
estnt[2]

betant = summary6ntdesc.mean[3:8]
estnt[1]


