module TestTrees

using Phylo
using DataFrames
using JuliaDB
using Query

using Compat.Test

species = ["Dog", "Cat", "Human"]
ntips = 10
df = DataFrame(species = species, count = [10, 20, 3])
observations = ["Dog", "Cat", "Dog", "Dog"]
jdb = table(@NT(species = observations, count = 1:4))
@testset "NamedTree()" begin
    ntn = NamedTree(ntips)
    @test length(nodefilter(isroot, ntn)) == ntips
    @test length(nodefilter(isleaf, ntn)) == ntips
    @test length(nodefilter(isinternal, ntn)) == 0
    @test validate(ntn)
    nt = NamedTree(species)
    @test_throws ErrorException noderoute(nt, "Dog", "Cat")
    @test_throws ErrorException branchroute(nt, "Dog", "Human")
    @test validate(nt)
    n = addnode!(nt)
    b1 = addbranch!(nt, n, "Dog", 2.0)
    b2 = addbranch!(nt, n, "Cat", 2.0)
    @test_throws ErrorException addbranch!(nt, n, "Human", 2.0)
    @test validate(nt)
    r = addnode!(nt)
    @test_throws ErrorException addbranch!(nt, r, "Potato", 2.0)
    b3 = addbranch!(nt, r, "Human", 5.0)
    b4 = addbranch!(nt, r, n, 3.0)
    @test maximum(distances(nt)) ≈ 10.0
    @test validate(nt)
    @test noderoute(nt, "Human", "Dog") == ["Human", r, n, "Dog"]
    @test branchroute(nt, "Human", "Dog") == [b3, b4, b1]
    nb = BinaryTree(DataFrame(name=["Human", "Cat"]))
    @test getleafinfo(nb) === getleafinfo(BinaryTree(nb))
    @test getleafinfo(nb) !== getleafinfo(BinaryTree(nb; copyinfo = true))
    @test_nowarn addnode!(nb, "Dog")
    @test_warn "LeafInfo names do not match actual leaves of tree" !validate(nb) || error("validate() should have returned false")
    np = PolytomousTree(DataFrame(name=["Human", "Cat"]))
    @test getleafinfo(np) === getleafinfo(PolytomousTree(np))
    @test getleafinfo(np) !== getleafinfo(PolytomousTree(np; copyinfo = true))
    @test_nowarn addnode!(np, "Dog")
    @test_warn "LeafInfo names do not match actual leaves of tree" !validate(np) || error("validate() should have returned false")
    @test_throws ErrorException PolytomousTree(np)
    @test getleafnames(nb) ⊆ getleafnames(np)
    @test getleafinfo(BinaryTree(df)) == df
    @test getleafinfo(PolytomousTree(df)) == df
end

@testset "BinaryTree()" begin
    btn = BinaryTree{DataFrame, Vector{Float64}}(ntips)
    @test length(nodefilter(isroot, btn)) == ntips
    @test length(nodefilter(isleaf, btn)) == ntips
    @test length(nodefilter(isinternal, btn)) == 0

    tmp = BinaryTree(DataFrame(names=species))
    @test validate(tmp)
    @test_nowarn addbranch!(tmp, addnode!(tmp), species[1])
    @test Set(getleafnames(tmp)) == Set(species)
    nt = BinaryTree{DataFrame, Vector{Float64}}(species)
    @test validate(nt)
    n = addnode!(nt)
    @test Set(getleafnames(nt)) == Set(species) ∪ Set([n])
    addbranch!(nt, n, "Dog", 2.0)
    addbranch!(nt, n, "Cat", 2.0)
    @test_throws ErrorException addbranch!(nt, n, "Human", 2.0)
    @test validate(nt)
    r = addnode!(nt)
    @test_throws ErrorException addbranch!(nt, r, "Potato", 2.0)
    addbranch!(nt, r, "Human", 4.0)
    addbranch!(nt, r, n, 2.0)
    @test validate(nt)
    setnoderecord!(nt, "Dog", [1.0])
    @test getnoderecord(nt, "Dog")[1] ≈ 1.0
    @test length(getnoderecord(nt, "Cat")) == 0
    nt2 = BinaryTree(nt)
    @test length(getnoderecord(nt2, "Dog")) == 0
    @test Set(getleafnames(nt2)) == Set(species)
    @test validate(nt)
    nt3 = BinaryTree(nt, empty=false)
    a=IOBuffer()
    b=IOBuffer()
    showall(a, nt)
    showall(b, nt3)
    @test a.data == b.data
    show(a, nt)
    show(b, nt3)
    @test a.data == b.data

    tdf = BinaryTree(df)
    @test nleaves(tdf) == 3
    tj = BinaryTree(jdb)
    @test nleaves(tj) == 2

    counts = @from line in getleafinfo(tj, "Dog") begin
        @select line.count
        @collect
    end
    @test sum(counts) == 8

    counts = @from line in getleafinfo(tj, "Cat") begin
        @select line.count
        @collect
    end
    @test sum(counts) == 2

    counts = @from line in getleafinfo(tdf, "Dog") begin
        @select line.count
        @collect
    end
    @test sum(counts) == 10

end

end
