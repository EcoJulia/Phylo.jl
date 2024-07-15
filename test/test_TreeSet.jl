# SPDX-License-Identifier: BSD-2-Clause

module TestTreeSet

using Phylo
using DataFrames

using Test

species = ["Dog", "Cat", "Human"]
ntips = 10
df = DataFrame(species = species, count = [10, 20, 3])
observations = ["Dog", "Cat", "Dog", "Dog"]
jdb = DataFrame(species = observations, count = [1, 2, 3, 4])

@testset "TreeSet" begin
    @test length(rand(Ultrametric(ntips), 10)) ==
          length(rand(Nonultrametric(ntips), 10)) == 10
    ts = rand(Ultrametric(ntips), 1:10)
    @test ts isa TreeSet
    @test isempty(getleafinfo(ts))
    @test length(gettrees(ts)) == length(gettreenames(ts)) == ntrees(ts) == 10
    @test length(getleafnames(ts)) == ntips
    for name in gettreenames(ts)
        @test name ∈ 1:10
        @test length(getbranchnames(ts[name])) == ntips * 2 - 2
        @test length(getbranchnames(ts)[name]) == ntips * 2 - 2
        @test length(getleafnames(ts[name])) == ntips
    end
    for tree in gettrees(ts)
        @test length(getnodenames(tree)) == ntips * 2 - 1
        @test length(getbranchnames(tree)) == ntips * 2 - 2
        @test nleaves(tree) == ntips
    end
    @test nleaves(ts) == ntips
    @test_throws ArgumentError getleafnames(TreeSet(Dict{String, NamedTree}(),
                                                    Dict{String,
                                                         Dict{String,
                                                              Any}}()))
    @test nleaves(rand(Ultrametric(df), 20)) == length(species)
    @test nleaves(rand(Ultrametric(jdb), 10)) ==
          length(unique(observations))
    @test getleafinfo(rand(Ultrametric(df))) == df
    @test getleafinfo(rand(Ultrametric(jdb), 2)) == jdb
    @test nleaves(rand(Nonultrametric(df), 20)) == length(species)
    @test nleaves(rand(Nonultrametric(jdb), 10)) ==
          length(unique(observations))
    @test getleafinfo(rand(Nonultrametric(df))) == df
    @test getleafinfo(rand(Nonultrametric(jdb), 2)) == jdb

    @test nodetype(typeof(ts)) == nodetype(typeof(first(gettrees(ts))))
    @test branchtype(typeof(ts)) == branchtype(typeof(first(gettrees(ts))))
end
end
