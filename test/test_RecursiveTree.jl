module TestRecursiveTree

using Phylo
using DataFrames

using Test
using IterableTables: getiterator

species = ["Dog", "Cat", "Human", "Potato"]
ntips = 10
df = DataFrame(species = species, count = [10, 20, 3, 31])
observations = ["Dog", "Cat", "Dog", "Dog"]
jdb = DataFrame(species = observations, count = 1:4)

@testset "RootedTree()" begin
    name = "internal"
    nts = RecursiveTree{OneRoot, String, Vector{Int}, Nothing, BinaryBranching, Float64, Nothing}(species)
    @test !renamenode!(nts, species[1], species[2])
    @test hasnode(nts, species[1])
    @test renamenode!(nts, species[1], "new " * species[1])
    @test_throws ErrorException renamenode!(nts, species[1], "new " * species[1])
    @test hasnode(nts, "new " * species[1])
    @test renamenode!(nts, "new " * species[1], species[1])
    @test Set(species) == Set(getleafnames(nts))

    rts = RootedTree(species)
    li = getleafinfo(rts)
    data = [1, 2]
    li[species[1]] = data
    @test !renamenode!(rts, species[1], species[2])
    @test renamenode!(rts, species[1], "new " * species[1])
    @test !haskey(li, species[1])
    @test haskey(li, "new " * species[1])
    @test li["new " * species[1]] == data
    @test_throws ErrorException renamenode!(rts, species[1], "new " * species[1])
    @test renamenode!(rts, "new " * species[1], species[1])
    @test haskey(li, species[1])
    @test li[species[1]] == data
    @test Set(species) == Set(getleafnames(rts))
    @test nodedatatype(typeof(rts)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(rts)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(rts)) ≡ Dict{String, Any}
    @test_nowarn createnode!(rts, name)
    @test createbranch!(rts, name, species[1]) ∈ getbranches(rts)
    
    rtdf = Phylo.ReT{OneRoot, DataFrame, BinaryBranching, Float64}(df)
    @test nodedatatype(typeof(rtdf)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(rtdf)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(rtdf)) ≡ DataFrame
    @test_nowarn createnode!(rtdf, name)
    @test createbranch!(rtdf, name, species[1]) ∈ getbranches(rtdf)
    @test createbranch!(rtdf, name, species[2]) ∈ getbranches(rtdf)
    @test_throws ErrorException createbranch!(rtdf, name, species[3])

    rtdfp = Phylo.ReT{OneRoot, DataFrame, PolytomousBranching, Float64}(df)
    @test nodedatatype(typeof(rtdfp)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(rtdfp)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(rtdfp)) ≡ DataFrame
    @test_nowarn createnode!(rtdfp, name)
    @test createbranch!(rtdfp, name, species[1]) ∈ getbranches(rtdfp)
    @test createbranch!(rtdfp, name, species[2]) ∈ getbranches(rtdfp)
    b = createbranch!(rtdfp, name, species[3])
    @test b ∈ getbranches(rtdfp)
    @test deletebranch!(rtdfp, b)
    @test createbranch!(rtdfp, name, species[3]) ∈ getbranches(rtdfp)
    @test_throws ErrorException createbranch!(rtdfp, name, species[2])
end

@testset "UnrootedTree()" begin
    name = "internal"
    urts = Phylo.ReTD{Unrooted, BinaryBranching, Float64}(species)
    @test nodedatatype(typeof(urts)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(urts)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(urts)) ≡ Dict{String, Any}
    @test_nowarn createnode!(urts, name)
    @test createbranch!(urts, name, species[1]) ∈ getbranches(urts)
    @test createbranch!(urts, name, species[2]) ∈ getbranches(urts)
    @test createbranch!(urts, name, species[3]) ∈ getbranches(urts)
    @test_throws ErrorException createbranch!(urts, name, species[4])

    urtsp = UnrootedTree(species)
    @test nodedatatype(typeof(urtsp)) ≡ Dict{String, Any}
    @test branchdatatype(typeof(urtsp)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(urtsp)) ≡ Dict{String, Any}
    @test_nowarn createnode!(urtsp, name)
    @test createbranch!(urtsp, name, species[1]) ∈ getbranches(urtsp)
    @test createbranch!(urtsp, name, species[2]) ∈ getbranches(urtsp)
    @test createbranch!(urtsp, name, species[3]) ∈ getbranches(urtsp)
    @test createbranch!(urtsp, name, species[4]) ∈ getbranches(urtsp)
    @test nroots(urtsp) == 0

    LB = RecursiveBranch{Unrooted, String, Vector{Int}, Nothing, BinaryBranching, Float64}
    LN = RecursiveNode{Unrooted, String, Vector{Int}, Nothing, BinaryBranching, Float64}
    rtjdb = RecursiveTree{Unrooted, String, Vector{Int}, Nothing, BinaryBranching, Float64, typeof(jdb)}(jdb)
    @test_throws MethodError renamenode!(rtjdb, observations[1], "new " * observations[1])
    @test nodedatatype(typeof(rtjdb)) ≡ Vector{Int}
    @test branchdatatype(typeof(rtjdb)) ≡ Nothing
    @test leafinfotype(typeof(rtjdb)) ≡ typeof(jdb)
    @test_nowarn createnode!(rtjdb, name)
    b = createbranch!(rtjdb, name, observations[1], data = nothing)
    @test b ∈ getbranches(rtjdb)
    @test deletebranch!(rtjdb, b)
    @test createbranch!(rtjdb, name, observations[1], data = nothing) ∈ getbranches(rtjdb)
    @test_nowarn Phylo.outputnode!(IOBuffer(), rtjdb, name, Phylo.CompactOutput(), Nothing)
    @test_nowarn Phylo.outputbranch!(IOBuffer(), rtjdb, first(getbranches(rtjdb)), Phylo.CompactOutput(), Nothing)
end

end
