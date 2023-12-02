module TestPlots
using Test

using Phylo
using Plots

@testset "Plots" begin
    tree = open(parsenewick, Phylo.path("hummingbirds.tree"))
    @test length(plot(tree).subplots) == 1
    trait = map_depthfirst((val, node) -> val + randn(), 0., tree, Float64)
    @test plot(tree, treetype = :fan, line_z = trait, linecolor = :RdYlBu, linewidth = 5, showtips = false).n == 1
    @test plot(tree,
        markersize = 10, markercolor = :steelblue, markerstrokecolor = :white,
        series_annotations = text.(1:nnodes(tree), 5, :center, :center, :white,
        tipfont = (4,))).init
end

end
