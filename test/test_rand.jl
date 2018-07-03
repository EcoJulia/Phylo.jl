module TestRand

using Phylo
using DataFrames

using Compat.Test

@testset "Nonultrametric()" begin
    # Create a 10 tip tree
    nu = Nonultrametric(10)
    @test eltype(nu) == NamedTree
    @test validate(rand(nu))
    @test Set(getleafnames(rand(nu))) == Set(getleafnames(rand(nu)))
    # Create a tree with named tips
    species = ["Dog", "Cat", "Human"]
    t = rand(Nonultrametric(species))
    @test ntrees(t) == 1
    @test validate(t)
    @test Set(getleafnames(t)) == Set(species)

    t2 = rand(Nonultrametric{BinaryTree{DataFrame, Vector{Float64}}}(species))
    @test length(getnoderecord(t2, species[1])) == 0
    t3 = rand(Nonultrametric{BinaryTree{DataFrame, Vector{String}}}(species))
    map(nodenameiter(t3)) do name
        setnoderecord!(t3, name, nodehistory(t3, name))
    end
    for name in nodenameiter(t3)
        @test all(getnoderecord(t3, name) .== nodehistory(t3, name))
    end
end

@testset "Ultrametric()" begin
    # Create a 10 tip tree
    u = Ultrametric(10)
    @test validate(rand(u))
    @test Set(getleafnames(rand(u))) == Set(getleafnames(rand(u)))
    tree = rand(u)
    @test ntrees(tree) == 1
    heights = map(x -> getheight(tree, x), getleafnames(tree))
    @test all(h -> h ≈ heights[1], heights)
    species = ["Dog", "Cat", "Human"]
    t = rand(Ultrametric(species))
    @test validate(t)
    @test Set(getleafnames(t)) == Set(species)

    numnodes = 50
    ul = Ultrametric{BinaryTree{DataFrame, Vector{Float64}}}(numnodes)
    u2 = rand(ul)
    @test eltype(ul) == BinaryTree{DataFrame, Vector{Float64}}
    @test length(nodefilter(isinternal, u2)) == numnodes - 2
end

@testset "TestSets" begin
    @test ntrees(rand(Ultrametric(10), 10)) == 10
    @test length(rand(Ultrametric(10), 10)) == 10
    names = ["One", "Two", "Three"]
    ts = rand(Nonultrametric(names), 20)
    @test length(ts) == 20
    @test ntrees(ts) == 20
    @test length(getleafnames(ts)) == length(names)
    @test length(getnodenames(ts)) == length(names) * 2 - 1
end
end
