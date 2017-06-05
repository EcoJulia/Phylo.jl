module TestTrees

using Phylo
using Base.Test

@testset "NamedTree()" begin
    nt = NamedTree(["Dog", "Cat", "Human"])
    @test validate(nt)
    n = addnode!(nt)
    @test !validate(nt)
    addbranch!(nt, n, "Dog", 2.0)
    addbranch!(nt, n, "Cat", 2.0)
    @test_throws ErrorException addbranch!(nt, n, "Human", 2.0)
    @test validate(nt)
    r = addnode!(nt)
    @test_throws ErrorException addbranch!(nt, r, "Potato", 2.0)
    @test !validate(nt)
    addbranch!(nt, r, "Human", 5.0)
    addbranch!(nt, r, n, 3.0)
    @test maximum(distance(nt)) ≈ 10.0
    @test validate(nt)
end

@testset "NodeTree()" begin
    nt = NodeTree(["Dog", "Cat", "Human"], nodetype=Vector{Float64})
    @test validate(nt)
    n = addnode!(nt)
    @test Set(getleafnames(nt)) == Set(["Dog", "Cat", "Human"])
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
    nt2 = NodeTree(nt)
    @test length(getnoderecord(nt2, "Dog")) == 0
    @test Set(getleafnames(nt2)) == Set(["Dog", "Cat", "Human"])
    setheight!(nt, "Dog", 4.0)
    setheight!(nt, "Cat", 5.0)
    setheight!(nt, "Human", 4.0)
    @test !validate(nt)
    setheight!(nt, "Cat", 4.0)
    @test validate(nt)
    setrootheight!(nt, 1.0)
    @test !validate(nt)
    setrootheight!(nt, 0.0)    
    @test validate(nt)
end

end
