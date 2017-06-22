module TestShow

using Phylo
using Base.Test
if !isdefined(Base.Test, Symbol("@test_nowarn"))
    # Ignore @test_nowarn unless it's there...
    macro test_nowarn(test)
    end
end

info("""
These tests only check that the show() and showall() commands
do not give warnings or errors, not that they are correct!
""")

ntips = 10
a = IOBuffer()
@testset "Nonultrametric()" begin
    nt = rand(Nonultrametric(ntips))
    @test_nowarn show(a, nt)
    @test_nowarn showall(a, nt)
    @test_nowarn show(a, first(NodeIterator(nt)))
    @test_nowarn show(a, first(BranchIterator(nt)))
    @test_nowarn show(a, first(NodeNameIterator(nt)) => first(NodeIterator(nt)))
    @test_nowarn show(a, first(BranchNameIterator(nt)) => first(NodeIterator(nt)))  
    bt = rand(Nonultrametric{BinaryTree{LeafInfo, Vector{String}}}(ntips))
    @test_nowarn show(a, bt)
    @test_nowarn showall(a, bt)
    @test_nowarn show(a, first(NodeIterator(bt)))
    @test_nowarn show(a, first(BranchIterator(bt)))
    @test_nowarn show(a, first(NodeNameIterator(bt)) => first(NodeIterator(bt)))
    @test_nowarn show(a, first(BranchNameIterator(bt)) => first(NodeIterator(bt)))
end
end
