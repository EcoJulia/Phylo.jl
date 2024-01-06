module TestLinkTree

using Phylo
using DataFrames

using Test
using IterableTables: getiterator

species = ["Dog", "Cat", "Human"]
ntips = 10
df = DataFrame(species = species, count = [10, 20, 3])
observations = ["Dog", "Cat", "Dog", "Dog"]
jdb = DataFrame(species = observations, count = 1:4)

@testset "RootedTree()" begin
    name = "internal"
    ltdf = Phylo.LT{OneRoot, DataFrame, Float64}(df)
    @test nodedatatype(typeof(ltdf)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(ltdf)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(ltdf)) ≡ DataFrame
    @test_nowarn createnode!(ltdf, name)
    @test createbranch!(ltdf, name, species[1]) ∈ getbranches(ltdf)
end

@testset "Unrooted Trees" begin
    name = "internal"
    LB = LinkBranch{Unrooted, String, Nothing, Float64}
    LN = LinkNode{Unrooted, String, Vector{Int}, LB}
    ltjdb = LinkTree{Unrooted, String, LN, LB, typeof(jdb)}(jdb)
    @test nodedatatype(typeof(ltjdb)) ≡ Vector{Int}
    @test branchdatatype(typeof(ltjdb)) ≡ Nothing
    @test leafinfotype(typeof(ltjdb)) ≡ typeof(jdb)
    @test_nowarn createnode!(ltjdb, name)
    @test createbranch!(ltjdb, name, observations[1], data = nothing) ∈
          getbranches(ltjdb)
    @test createbranch!(ltjdb, getnode(ltjdb, name),
                        getnode(ltjdb, observations[2]), data = nothing) ∈
          getbranches(ltjdb)
end

end
