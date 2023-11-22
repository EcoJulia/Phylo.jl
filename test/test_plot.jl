module TestPlots
using Test

using Phylo
using Plots

@testset "Plots" begin
    tree = open(parsenewick, Phylo.path("hummingbirds.tree"))
    @test all(getproperty.([plot(tree), plot(tree, treetype = :fan)], :n) .== 1)
end

end
