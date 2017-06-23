module TestRand

using Phylo
using Base.Test

@testset "Nonultrametric()" begin
    # Create a 10 tip tree
    nu = Nonultrametric(10)
    @test eltype(nu) == NamedTree
    @test validate(rand(nu))
    @test Set(getleafnames(rand(nu))) == Set(getleafnames(rand(nu)))
    # Create a tree with named tips
    species = ["Dog", "Cat", "Human"]
    t = rand(Nonultrametric(species))
    @test validate(t)
    @test Set(getleafnames(t)) == Set(species)

    t2 = rand(Nonultrametric{BinaryTree{LeafInfo, Vector{Float64}}}(species))
    @test length(getnoderecord(t2, species[1])) == 0
    t3 = rand(Nonultrametric{BinaryTree{LeafInfo, Vector{String}}}(species))
    map(node -> setnoderecord!(t3, node, nodehistory(t3, node)), NodeNameIterator(t3))
    for name in NodeNameIterator(t3)
        @test all(getnoderecord(t3, name) .== nodehistory(t3, name))
    end
end

@testset "Ultrametric()" begin
    # Create a 10 tip tree
    u = Ultrametric(10)
    @test validate(rand(u))
    @test Set(getleafnames(rand(u))) == Set(getleafnames(rand(u)))
    tree = rand(u)
    heights = map(x -> getheight(tree, x), getleafnames(tree))
    @test all(h -> h ≈ heights[1], heights)
    species = ["Dog", "Cat", "Human"]
    t = rand(Ultrametric(species))
    @test validate(t)
    @test Set(getleafnames(t)) == Set(species)

    numnodes = 50
    ul = Ultrametric{BinaryTree{LeafInfo, Vector{Float64}}}(numnodes)
    u2 = rand(ul)
    @test eltype(ul) == BinaryTree{LeafInfo, Vector{Float64}}
    @test length(NodeIterator(isinternal, u2)) == numnodes - 2
end
end
