module TestRand

using Phylo
using Base.Test

@testset "Nonultrametric()" begin
    # Create a 10 tip tree
    nu = Nonultrametric(10)
    @test validate(rand(nu))
    @test Set(getleafnames(rand(nu))) == Set(getleafnames(rand(nu)))
    # Create a tree with named tips
    species = ["Dog", "Cat", "Human"]
    t = rand(Nonultrametric(species))
    @test validate(t)
    @test Set(getleafnames(t)) == Set(species)
end

@testset "Ultrametric()" begin
    # Create a 10 tip tree
    u = Ultrametric(10)
    @test validate(rand(u))
    @test Set(getleafnames(rand(u))) == Set(getleafnames(rand(u)))
    tree = rand(u)
    heights = map(x -> getheight(tree, x), getleafnames(tree))
    @test all(h -> h ≈ heights[1], heights)
end
end
