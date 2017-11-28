module TestTrees

using Phylo
using Compat.Test

species = ["Dog", "Cat", "Human"]
ntips = 10
@testset "NamedTree()" begin
    ntn = NamedTree(ntips)
    @test length(nodefilter(isroot, ntn)) == ntips
    @test length(nodefilter(isleaf, ntn)) == ntips
    @test length(nodefilter(isinternal, ntn)) == 0
    @test validate(ntn)
    nt = NamedTree(species)
    @test isnull(noderoute(nt, "Dog", "Cat"))
    @test isnull(branchroute(nt, "Dog", "Human"))
    @test validate(nt)
    n = addnode!(nt)
    @test_warn "Leaf names do not match actual leaves of tree" !validate(nt) || error("validate() should have returned false")
    b1 = addbranch!(nt, n, "Dog", 2.0)
    b2 = addbranch!(nt, n, "Cat", 2.0)
    @test_throws ErrorException addbranch!(nt, n, "Human", 2.0)
    @test validate(nt)
    r = addnode!(nt)
    @test_throws ErrorException addbranch!(nt, r, "Potato", 2.0)
    @test_warn "Leaf names do not match actual leaves of tree" !validate(nt) || error("validate() should have returned false")
    b3 = addbranch!(nt, r, "Human", 5.0)
    b4 = addbranch!(nt, r, n, 3.0)
    @test maximum(distances(nt)) ≈ 10.0
    @test validate(nt)
    @test get(noderoute(nt, "Human", "Dog")) == ["Human", r, n, "Dog"]
    @test get(branchroute(nt, "Human", "Dog")) == [b3, b4, b1]
end

@testset "BinaryTree()" begin
    btn = BinaryTree{LeafInfo, Vector{Float64}}(ntips)
    @test length(nodefilter(isroot, btn)) == ntips
    @test length(nodefilter(isleaf, btn)) == ntips
    @test length(nodefilter(isinternal, btn)) == 0

    nt = BinaryTree{LeafInfo, Vector{Float64}}(species)
    @test validate(nt)
    n = addnode!(nt)
    @test Set(getleafnames(nt)) == Set(species)
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
    setheight!(nt, "Dog", 4.0)
    setheight!(nt, "Cat", 5.0)
    setheight!(nt, "Human", 4.0)
    @test_warn "does not match branches" !validate(nt) || error("validate() should have returned false")
    setheight!(nt, "Cat", 4.0)
    @test validate(nt)
    setrootheight!(nt, 1.0)
    @test_warn "does not match branches" !validate(nt) || error("validate() should have returned false")
    setrootheight!(nt, 0.0)    
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
end

end
