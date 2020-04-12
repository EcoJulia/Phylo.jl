module TestLinkTree

using Phylo
using DataFrames
using JuliaDB

using Test
using IterableTables: getiterator

species = ["Dog", "Cat", "Human"]
ntips = 10
df = DataFrame(species = species, count = [10, 20, 3])
observations = ["Dog", "Cat", "Dog", "Dog"]
jdb = table((species = observations, count = 1:4))

@testset "LinkTree()" begin
    lts = RootedTree(species)
    @test nodedatatype(typeof(lts)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(lts)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(lts)) ≡ Dict{String, Any}
    
    ltdf = Phylo.LT{OneRoot, DataFrame, Float64}(df)
    @test nodedatatype(typeof(ltdf)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(ltdf)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(ltdf)) ≡ DataFrame

    rt = Unrooted
    lb = LinkBranch{rt, String, Nothing}
    ln = LinkNode{rt, String, Vector{Int}, lb} 
    ltjdb = LinkTree{rt, String, ln, lb, typeof(jdb)}(jdb)
    @test nodedatatype(typeof(ltjdb)) ≡ Vector{Int}
    @test branchdatatype(typeof(ltjdb)) ≡ Nothing
    @test leafinfotype(typeof(ltjdb)) ≡ typeof(jdb)
end

end
