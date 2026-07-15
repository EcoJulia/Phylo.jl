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

# Bayesian inference on Myrtaceae data, all traits 

# Load the trees
include("read_tree_13.jl");

@model function βσ_multthreepoint(tree, z) 

    β ~ MvNormal(zeros(13), 1.0 * I)
    #σ ~ InverseWishart(13.0, PDMat(I(13)*1.0))
    
    σ_scale ~ filldist(truncated(TDist(3); lower = 0), 13)  # positive half of t dist
    #Lcorr ~ LKJCholesky(13, 1.0)
    σ_corr ~ LKJ(13,1)
    #σ = diagm(σ_scale) * (Lcorr.L * Lcorr.L') * diagm(σ_scale)
    σ = diagm(σ_scale) * σ_corr * diagm(σ_scale)

    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end


n_samples = 100_000
n_warmup = 10_000

model = βσ_multthreepoint(tree, z);
#init = InitFromParams((β = zeros(13), σ_scale = ones(13), Lcorr = rand(LKJCholesky(13, 0.1))))
init = InitFromParams((β = zeros(13), σ_scale = ones(13), Lcorr = Matrix(I*1.0, 13, 13)))

chn = sample(model, NUTS(; adtype=AutoReverseDiff(; compile=true)), n_samples; 
                initial_params = init,
                num_warmup = n_warmup)
# LKJCholesky ~40 hours
# LKJ ~90 hours
println(describe(chn))

summary = describe(chn)
summary_stats = summary[1]
summarydf = DataFrame(summary_stats)

plot(chn6_1["β[1]"])
plot(chn6_1["σ_scale[1]"])
plot(chn6_1["Lcorr.L[6, 1]"])

est = estimaterates(tree6, trait)

sig_scale = summarydf.mean[7:12]

v = summarydf.mean[13:33]
corr = LowerTriangular(reshape([v; zeros(6^2 - length(v))], 6, 6))

corrl = zeros(6, 6)
corrl[tril!(trues(6, 6))] = v

sigma = diagm(sig_scale) * (corrl * corrl') * diagm(sig_scale)
est[2]

beta = summarydf.mean[1:6]
est[1]