module TestInference_wrapped

using Test
using Phylo
using RCall
using DataFrames

global skipR = !rcopy(R"require(ape)")
global skipR = !rcopy(R"require(phylolm)")

@testset "Compare estimaterates output to phylolm" begin
    
    R"
    library('ape')
    library('phylolm')"

    jtree = open(parsenewick, Phylo.path("hummingbirds.tree"))
    R"rtree = ape::read.tree(file = 'hummingbirds.tree')" #I dont think this is the way I should load in the tree to R. Tried like in other files but I couldnt get it to load into phylolm.

    # Create dataframe with leafnames and random trait values
    species = getleafnames(jtree)
    data = 1000 .* rand(length(species))
    dat = DataFrame(species = species, data = data)
    @rput dat

    # Save data on leaves so can test on estimaterates 2
    setnodedata!.(jtree, species, "trait", data)

    # Species need to be the row names in R
    R"row.names(dat) <- dat$species"

    rfit = rcopy(R"phylolm(data~1,data=dat,phy=rtree)")
    jfit = estimaterates(jtree, "trait")

    species = getleaves(jtree)
    data = 1000 .* rand(length(species))

    @test rfit[:coefficients] ≈ jfit[1]
    @test rfit[:sigma2] ≈ jfit[2]
    @test rfit[:logLik] ≈ -jfit[3]
end

end