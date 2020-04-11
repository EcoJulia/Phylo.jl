module TestRand

using Phylo
using DataFrames
using Random

using Test

@testset "Nonultrametric()" begin
    # Create a 10 tip tree
    nu = Nonultrametric(10)
    @test eltype(nu) == RootedTree
    @test validate!(rand(nu))
    @test Set(getleafnames(rand(nu))) == Set(getleafnames(rand(nu)))
    # Create a tree with named tips
    species = ["Dog", "Cat", "Human"]
    t = rand(Nonultrametric(species))
    @test ntrees(t) == 1
    @test validate!(t)
    @test Set(getleafnames(t)) == Set(species)

    t2 = rand(Nonultrametric{BinaryTree{OneRoot, DataFrame, Vector{Float64}}}(species))
    @test length(getnodedata(t2, species[1])) == 0
    t3 = rand(Nonultrametric{BinaryTree{OneRoot, DataFrame, Vector{String}}}(species))
    map(nodenameiter(t3)) do name
        setnodedata!(t3, name, nodehistory(t3, name))
    end
    for name in nodenameiter(t3)
        @test all(getnodedata(t3, name) .== nodehistory(t3, name))
    end
    d = rand(BrownianTrait(t, "trait", σ² = 1.0))
    @test d ≢ t
    @test species ⊆ collect(keys(d))
    @test eltype(values(d)) ≡ typeof(1.0)
    t4 = rand!(BrownianTrait(t, "trait64", σ² = 1.0), t)
    @test t4 ≡ t
    @test typeof(getnodedata(t4, species[1])["trait64"]) ≡ typeof(1.0)
    t5 = rand!(BrownianTrait(t, "trait32", 0.0f0, σ² = 1.0f0), t)
    @test typeof(getnodedata(t5, species[1])["trait32"]) ≡ typeof(1.0f0)
    @test t5 ≡ t
end

@testset "Ultrametric()" begin
    # Create a 10 tip tree
    u = Ultrametric(10)
    @test validate!(rand(u))
    @test Set(getleafnames(rand(u))) == Set(getleafnames(rand(u)))
    tree = rand(u)
    @test ntrees(tree) == 1
    heights = map(x -> getheight(tree, x), getleafnames(tree))
    @test all(h -> h ≈ heights[1], heights)
    species = ["Dog", "Cat", "Human"]
    t = rand(Ultrametric(species))
    @test validate!(t)
    @test Set(getleafnames(t)) == Set(species)

    numnodes = 50
    ul = Ultrametric{BinaryTree{OneRoot, DataFrame, Vector{Float64}}}(numnodes)
    u2 = rand(ul)
    @test eltype(ul) == BinaryTree{OneRoot, DataFrame, Vector{Float64}}
    @test length(nodefilter(isinternal, u2)) == numnodes - 2

    d = rand(BrownianTrait(t, "trait", σ² = 1.0))
    @test species ⊆ collect(keys(d))
    @test eltype(values(d)) ≡ typeof(1.0)
    t4 = rand!(BrownianTrait(t, "trait64", σ² = 1.0), t)
    @test t4 ≡ t
    @test typeof(getnodedata(t4, species[1])["trait64"]) ≡ typeof(1.0)
    t5 = rand!(BrownianTrait(t, "trait32", 0.0f0, σ² = 1.0f0), t)
    @test typeof(getnodedata(t5, species[1])["trait32"]) ≡ typeof(1.0f0)
    @test t5 ≡ t
end

@testset "TestSets" begin
    @test ntrees(rand(Ultrametric(10), 10)) == 10
    @test length(rand(Ultrametric(10), 10)) == 10
    names = ["One", "Two", "Three"]
    ts = rand(Nonultrametric(names), 20)
    @test length(ts) == 20
    @test ntrees(ts) == 20
    @test length(getleafnames(ts)) == length(names)
    @test length(getnodenames(first(gettrees(ts)))) .== length(names) * 2 - 1
end
end
