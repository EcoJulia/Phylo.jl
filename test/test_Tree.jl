module TestTrees

using Phylo
using DataFrames
using JuliaDB

using Compat.Test
using IterableTables: getiterator

species = ["Dog", "Cat", "Human"]
ntips = 10
df = DataFrame(species = species, count = [10, 20, 3])
observations = ["Dog", "Cat", "Dog", "Dog"]
@static if VERSION < v"0.7.0-"
    jdb = table(@NT(species = observations, count = 1:4))
else
    jdb = table((species = observations, count = 1:4))
end
@testset "NamedBinaryTree()" begin
    ntn = NamedBinaryTree(ntips)
    @test length(nodefilter(isroot, ntn)) == ntips
    @test length(nodefilter(isleaf, ntn)) == ntips
    @test length(nodefilter(isinternal, ntn)) == 0
    @test validate!(ntn)
    nt = NamedBinaryTree(species)
    @test_throws Exception noderoute(nt, "Dog", "Cat")
    @test_throws Exception branchroute(nt, "Dog", "Human")
    @test validate!(nt)
    n = getnodename(nt, createnode!(nt))
    b1 = createbranch!(nt, n, "Dog", 2.0)
    b2 = createbranch!(nt, n, "Cat", 2.0)
    @test_throws Exception createbranch!(nt, n, "Human", 2.0)
    @test validate!(nt)
    r = getnodename(nt, createnode!(nt))
    @test_throws Exception createbranch!(nt, r, "Potato", 2.0)
    b3 = createbranch!(nt, r, "Human", 5.0)
    b4 = createbranch!(nt, r, n, 3.0)
    @test maximum(distances(nt)) ≈ 10.0
    @test validate!(nt)
    @test noderoute(nt, "Human", "Dog") == ["Human", r, n, "Dog"]
    @test branchroute(nt, "Human", "Dog") == [b3, b4, b1]
    nb = BinaryTree{OneRoot}(DataFrame(name=["Human", "Cat"]))
    @test_throws Exception BinaryTree(nb)
    nb = BinaryTree(DataFrame(name=["Human", "Cat"]))
    @test_throws Exception BinaryTree(nb)
    nb = BinaryTree(DataFrame(name=["Human"]))
    @test_nowarn BinaryTree(nb)
    nb = BinaryTree{ManyRoots}(DataFrame(name=["Human", "Cat"]))
    @test getleafinfo(nb) === getleafinfo(BinaryTree(nb))
    @test getleafinfo(nb) !== getleafinfo(BinaryTree(nb; copyinfo = true))
    @test_nowarn createnode!(nb, "Dog")
    if VERSION < v"0.7.0-"
        @test_warn "LeafInfo names do not match actual leaves of tree" !validate!(nb) || error("validate!() should have returned false")
    else
        @test !validate!(nb)
    end
    np = PolytomousTree{OneRoot}(DataFrame(name=["Human", "Cat"]))
    @test_throws Exception PolytomousTree(np)
    np = PolytomousTree(DataFrame(name=["Human", "Cat"]))
    @test_throws Exception BinaryTree(np)
    np = PolytomousTree(DataFrame(name=["Human"]))
    @test_nowarn PolytomousTree(np)
    np = PolytomousTree{ManyRoots}(DataFrame(name=["Human", "Cat"]))
    @test getleafinfo(np) !== getleafinfo(PolytomousTree(np; copyinfo = true))
    @test_nowarn createnode!(np, "Dog")
    if VERSION < v"0.7.0-"
        @test_warn "LeafInfo names do not match actual leaves of tree" !validate!(np) || error("validate!() should have returned false")
    else
        @test !validate!(np)
    end
    @test_throws Exception PolytomousTree(np)
    @test getleafnames(nb) ⊆ getleafnames(np)
    @test getleafinfo(BinaryTree(df)) == df
    @test getleafinfo(PolytomousTree(df)) == df
    @test branchdatatype(typeof(np)) ≡ Nothing
    @test nodedatatype(typeof(np)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(np)) ≡ DataFrame
end

@testset "BinaryTree()" begin
    btn = BinaryTree{OneRoot, DataFrame, Vector{Float64}}(ntips)
    @test length(nodefilter(isroot, btn)) == ntips
    @test length(nodefilter(isleaf, btn)) == ntips
    @test length(nodefilter(isinternal, btn)) == 0

    tmp = BinaryTree{ManyRoots}(DataFrame(names=species))
    @test validate!(tmp)
    tmp = BinaryTree(DataFrame(names=species))
    @test !validate!(tmp)
    names = createnodes!(tmp, 2)
    @test_nowarn createbranch!(tmp, names[1], species[1])
    @test_nowarn createbranch!(tmp, names[1], species[2])
    @test_nowarn createbranch!(tmp, names[2], species[3])
    @test_nowarn createbranch!(tmp, names[2], names[1])
    @test Set(getleafnames(tmp)) == Set(species)
    @test validate!(tmp)
    nt = BinaryTree{ManyRoots, DataFrame, Vector{Float64}}(species)
    @test validate!(nt)
    n = getnodename(nt, createnode!(nt))
    @test Set(getleafnames(nt)) == Set(species) ∪ Set([n])
    createbranch!(nt, n, "Dog", 2.0)
    createbranch!(nt, n, "Cat", 2.0)
    @test_throws Exception createbranch!(nt, n, "Human", 2.0)
    @test validate!(nt)
    r = getnodename(nt, createnode!(nt))
    @test_throws Exception createbranch!(nt, r, "Potato", 2.0)
    createbranch!(nt, r, "Human", 4.0)
    createbranch!(nt, r, n, 2.0)
    @test validate!(nt)
    setnodedata!(nt, "Dog", [1.0])
    @test getnodedata(nt, "Dog")[1] ≈ 1.0
    @test length(getnodedata(nt, "Cat")) == 0
    nt2 = BinaryTree(nt)
    @test length(getnodedata(nt2, "Dog")) == 0
    @test Set(getleafnames(nt2)) == Set(species)
    @test validate!(nt)
    nt3 = BinaryTree(nt)
    nt4 = BinaryTree(nt3)
    a=IOBuffer()
    b=IOBuffer()
    show(a, nt3)
    show(b, nt4)
    @test all(a.data .== b.data)

    tdf = BinaryTree(df)
    @test nleaves(tdf) == 3
    tj = BinaryTree(jdb)
    @test nleaves(tj) == 2

    lij = getleafinfo(tj, "Dog")
    @test sum(map(line -> line.count, lij)) == 8

    lij = getleafinfo(tj, "Cat")
    @test sum(map(line -> line.count, lij)) == 2

    lidf = getleafinfo(tdf, "Dog")
    @test sum(map(line -> line.count, lidf)) == 10

    @test branchdatatype(typeof(tdf)) ≡ Nothing
    @test nodedatatype(typeof(tdf)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(tdf)) ≡ DataFrame
    @test branchdatatype(typeof(tj)) ≡ Nothing
    @test nodedatatype(typeof(tj)) ≡ Dict{String, Any}
    @test leafinfotype(typeof(tj)) <: IndexedTable
end

end
