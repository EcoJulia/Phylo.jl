module TestInference_wrapped

using Test
using Phylo
using RCall
using DataFrames

global skipR = !(rcopy(R"require(ape)") && rcopy(R"require(phylolm)"))

@testset "Compare estimaterates output to phylolm" begin
    jtree = nothing
    @test_nowarn jtree = open(parsenewick, Phylo.path("hummingbirds.tree"))

    # Create dataframe with leafnames and random trait values
    species = getleafnames(jtree)
    data = 1000 .* rand(length(species))
    dat = DataFrame(species = species, data = data)

    # Save data on leaves so can test on estimaterates 2
    setnodedata!.(jtree, species, "trait", data)

    jfit = nothing
    @test_nowarn jfit = estimaterates(jtree, "trait")

    species = getleaves(jtree)
    data = 1000 .* rand(length(species))

    if !skipR
        rtree = rcall(Symbol("read.tree"), Phylo.path("hummingbirds.tree"))
        @rput rtree
        @rput dat

        # Species need to be the row names in R
        R"row.names(dat) <- dat$species"

        rfit = rcopy(R"phylolm(data ~ 1, data = dat, phy = rtree)")

        @test rfit[:coefficients] ≈ jfit[1]
        @test rfit[:sigma2] ≈ jfit[2]
        @test rfit[:logLik] ≈ -jfit[3]
    end
end

end
