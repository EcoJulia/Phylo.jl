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
    lts = RootedTree(species)
    @test nodedatatype(typeof(lts)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(lts)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(lts)) ≡ Dict{String, Any}
    @test_nowarn createnode!(lts, name)
    @test createbranch!(lts, name, species[1]) ∈ getbranches(lts)
    
    ltdf = Phylo.LT{OneRoot, DataFrame, Float64}(df)
    @test nodedatatype(typeof(ltdf)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(ltdf)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(ltdf)) ≡ DataFrame
    @test_nowarn createnode!(ltdf, name)
    @test createbranch!(ltdf, name, species[1]) ∈ getbranches(ltdf)
end

@testset "UnrootedTree()" begin
    name = "internal"
    urts = UnrootedTree(species)
    @test nodedatatype(typeof(urts)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(urts)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(urts)) ≡ Dict{String, Any}
    @test_nowarn createnode!(urts, name)
    @test createbranch!(urts, name, species[1]) ∈ getbranches(urts)

    RT = Unrooted
    LB = LinkBranch{RT, String, Nothing, Float64}
    LN = LinkNode{RT, String, Vector{Int}, LB}
    ltjdb = LinkTree{RT, String, LN, LB, typeof(jdb)}(jdb)
    @test nodedatatype(typeof(ltjdb)) ≡ Vector{Int}
    @test branchdatatype(typeof(ltjdb)) ≡ Nothing
    @test leafinfotype(typeof(ltjdb)) ≡ typeof(jdb)
    @test_nowarn createnode!(ltjdb, name)
    @test createbranch!(ltjdb, name, observations[1], data = nothing) ∈ getbranches(ltjdb)
end

end
