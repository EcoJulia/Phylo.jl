module TestPlots
using Test

using Phylo
using Plots
using Random

@testset "Plots" begin
    tree = open(parsenewick, Phylo.path("hummingbirds.tree"));
    @test length(plot(tree).subplots) == 1
    trait = map_depthfirst((val, node) -> val + randn(), 0., tree, Float64)
    @test plot(sort!(tree), treetype = :fan, line_z = trait, linecolor = :RdYlBu, linewidth = 5, showtips = false).n == 1
    @test plot(tree, markersize = 10, markercolor = :steelblue,
               markerstrokecolor = :white,
               series_annotations = text.(1:nnodes(tree), 5, :center, :center, :white),
                                          tipfont = (4,)).init

    @enum TemperatureTrait lowTempPref midTempPref highTempPref
    tempsampler = SymmetricDiscreteTrait(tree, TemperatureTrait, 0.4, "Temperature")
    rand!(tempsampler, tree)
        
    ## and plot it
    @test plot(tree, showtips = true,
               marker_group = "Temperature",
               legend = :topleft, msc = :white, treetype = :fan,
               c = [:red :blue :green]).init
end

end
