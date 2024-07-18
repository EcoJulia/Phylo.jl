# SPDX-License-Identifier: BSD-2-Clause

module TestInference_wrapped

using Test
using Phylo
using RCall
using DataFrames
using DynamicPPL
using ForwardDiff

global skipR = !(rcopy(R"require(ape)") && rcopy(R"require(phylolm)"))

@testset "Compare estimaterates output to phylolm" begin
    jtree = open(f -> parsenewick(f, TraitTree{1}),
                 Phylo.path("hummingbirds.tree"))

    # Create dataframe with leafnames and random trait values
    species = getleafnames(jtree)
    data = 1000 .* rand(length(species))
    dat = DataFrame(species = species, data = data)

    # Save data on leaves so can test on estimaterates 2
    traits = Phylo.traitdata.(Union{Float64,
                                    ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag,
                                                                     Float64},
                                                     Float64, 3}}, Ref("trait"),
                              dat.data)
    setnodedata!.(jtree, dat.species, traits)

    jfit = nothing
    @test_nowarn jfit = estimaterates(jtree, ["trait"])

    if !skipR
        rtree = rcall(Symbol("read.tree"), Phylo.path("hummingbirds.tree"))
        @rput rtree
        @rput dat

        # Species need to be the row names in R
        R"row.names(dat) <- dat$species"

        rfit = rcopy(R"phylolm(data ~ 1, data = dat, phy = rtree)")

        @test rfit[:coefficients] ≈ jfit[1][1]
        @test rfit[:sigma2] ≈ jfit[2][1]
        @test rfit[:logLik] ≈ -jfit[3][1]
    end
end

end
