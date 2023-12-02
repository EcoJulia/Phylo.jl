module TestPlots
using Test

using Phylo
using Plots

@testset "Plots" begin
    tree = open(parsenewick, Phylo.path("hummingbirds.tree"))
    @test plot(tree).n == 1
    trait = map_depthfirst((val, node) -> val + randn(), 0., tree, Float64)
    @test plot(tree, treetype = :fan, line_z = trait, linecolor = :RdYlBu, linewidth = 5, showtips = false).n == 1
end

end
