# SPDX-License-Identifier: BSD-2-Clause

using Turing
using Phylo
using Distributions
using Random
using BenchmarkTools
using DataFrames
using LinearAlgebra
using Plots
using CSV
using ForwardDiff
using DynamicPPL
using PDMats


# How code run time scales with number of traits

@model function βσ()
    β ~ Uniform(-100, 1000)
    σ ~ Uniform(0, 100)
    return β, σ
end

@model function βσ_threepoint(tree, z) 
    @submodel β, σ = βσ()
    z ~ Phylo.MyDist2(σ, β, tree) 
    return nothing
end

@model function βσλ_threepoint(tree, z, upper = 1.0) 
    @submodel β, σ = βσ()
    λ ~ Uniform(0, upper)
    z ~ Phylo.MyDist3(σ, β, λ, tree)
    return nothing
end

@model function βσ_covariance(z, C)
    @submodel β, σ = βσ()
    z ~ MvNormal(β * ones(length(z)), σ * C)
    return nothing
end

@model function βσ_mult2()

    β_1 ~ Uniform(-100, 1000)
    β_2 ~ Uniform(-100, 1000)
    σ_1 ~ Uniform(0, 100)

    β = [β_1, β_2]
   
    σ ~ InverseWishart(2.0, PDMat(I(2)*1.0))

    return β, σ
end


@model function βσ_multthreepoint2(tree, z) 
    @submodel β, σ = βσ_mult2()
    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end

@model function βσ_covariance2(z, C)
    @submodel β, σ = βσ_mult2()
    n = dim(C)
    z ~ MvNormal(repeat(β, outer=n), kron(σ, C))
    return nothing
end

@model function βσ_mult3()
    β_1 ~ Uniform(-100, 1000)
    β_2 ~ Uniform(-100, 1000)
    β_3 ~ Uniform(-100, 1000)

    σ ~ InverseWishart(3.0, PDMat(I(3)*1.0))

    β = [β_1, β_2, β_3]

    return β, σ
end


@model function βσ_multthreepoint3(tree, z) 
    @submodel β, σ = βσ_mult3()
    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end

@model function βσ_mult4()
    β_1 ~ Uniform(-100, 1000)
    β_2 ~ Uniform(-100, 1000)
    β_3 ~ Uniform(-100, 1000)
    β_4 ~ Uniform(-100, 1000)

    β = [β_1, β_2, β_3, β_4]

    σ ~ InverseWishart(4.0, PDMat(I(4)*1.0))

    return β, σ
end

@model function βσ_multthreepoint4(tree, z) # z needs to be for leaves in postorder
    @submodel β, σ = βσ_mult4()
    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end

@model function βσ_mult5()
    β_1 ~ Uniform(-100, 1000)
    β_2 ~ Uniform(-100, 1000)
    β_3 ~ Uniform(-100, 1000)
    β_4 ~ Uniform(-100, 1000)
    β_5 ~ Uniform(-100, 1000)

    β = [β_1, β_2, β_3, β_4, β_5]

    σ ~ InverseWishart(5.0, PDMat(I(5)*1.0))

    return β, σ
end

@model function βσ_multthreepoint5(tree, z) # z needs to be for leaves in postorder
    @submodel β, σ = βσ_mult5()
    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end

@model function βσ_mult13()
    β_1 ~ Uniform(-100, 1000)
    β_2 ~ Uniform(-100, 1000)
    β_3 ~ Uniform(-100, 1000)
    β_4 ~ Uniform(-100, 1000)
    β_5 ~ Uniform(-100, 1000)
    β_6 ~ Uniform(-100, 1000)
    β_7 ~ Uniform(-100, 1000)
    β_8 ~ Uniform(-100, 1000)
    β_9 ~ Uniform(-100, 1000)
    β_10 ~ Uniform(-100, 1000)
    β_11 ~ Uniform(-100, 1000)
    β_12 ~ Uniform(-100, 1.0e8)
    β_13 ~ Uniform(-100, 1000)

    β = [β_1, β_2, β_3, β_4, β_5, β_6, β_7, β_8, β_9, β_10, β_11, β_12, β_13]

    σ ~ InverseWishart(13.0, PDMat(I(13)*1.0))

    return β, σ
end

@model function βσ_multthreepoint13(tree, z) # z needs to be for leaves in postorder
    @submodel β, σ = βσ_mult13()
    z ~ Phylo.MyDist4(σ, β, tree) 
    return nothing
end

# Load the trees
include("read_trees.jl");


## SCALING FOR MULTIPLE TRAITS 
n_samples = 1_000

# 1 trait
model_tp1 = βσ_threepoint(bigtree, z1);
spltp1 = sample(model_tp1, HMC(0.01, 5), n_samples)#; initial_params = [290, 8])
# old 21.83 seconds

# 1 trait w/ lambda
model_tpl = βσλ_threepoint(bigtree, dat.tmin);
spltpl = sample(model_tpl, HMC(0.01, 5), n_samples)
# old 39.85 seconds

model_c1 = βσ_covariance(z1, C);
splc1 = sample(model_c1, HMC(0.01, 5), n_samples)
# old 540.84 seconds

# 2 traits
model_tp2 = βσ_multthreepoint2(bigtree2, z2);
spl2 = sample(model_tp2, HMC(0.01, 5), n_samples)
# old 35.4 seconds

model_c2 = βσ_covariance2(z2, C);
splc2 = sample(model_c2, HMC(0.01, 5), n_samples)
# old estimated 2 hours

# 3 traits
model_tp3 = βσ_multthreepoint3(bigtree3, z3);
spl3 = sample(model_tp3, HMC(0.01, 5), n_samples)
# old 34.6 seconds


# 4 traits
model_tp4 = βσ_multthreepoint4(bigtree4, z4);
spl4 = sample(model_tp4, HMC(0.01, 5), n_samples)
# old 50.06 seconds

# 5 traits
model_tp5 = βσ_multthreepoint5(bigtree5, z5);
spl5 = sample(model_tp5, HMC(0.01, 5), n_samples)
# old 65.64 seconds

# 13 traits
model_tp13 = βσ_multthreepoint13(bigtree13, z13);
spl13 = sample(model_tp13, HMC(0.01, 5), n_samples)
# PosDefException: matrix is not positive definite; Factorization failed.

