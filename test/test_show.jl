module TestShow

using Phylo
using DataFrames

using Compat.Test
using Compat: @info

@info """
These tests only check that the show() commands
do not give warnings or errors, not that they are correct!
"""

ntips = 10
a = IOBuffer()
@testset "Nonultrametric()" begin
    nt = rand(Nonultrametric(ntips))
    @test_nowarn show(a, nt)
    @test_nowarn show(a, first(nodeiter(nt)))
    @test_nowarn show(a, first(branchiter(nt)))
    @test_nowarn show(a, first(nodenameiter(nt)) => first(nodeiter(nt)))
    @test_nowarn show(a, first(branchnameiter(nt)) => first(nodeiter(nt)))
    bt = rand(Nonultrametric{BinaryTree{DataFrame, Vector{String}}}(ntips))
    @test_nowarn show(a, bt)
    @test_nowarn show(a, first(nodeiter(bt)))
    @test_nowarn show(a, first(branchiter(bt)))
    @test_nowarn show(a, first(nodenameiter(bt)) => first(nodeiter(bt)))
    @test_nowarn show(a, first(branchnameiter(bt)) => first(nodeiter(bt)))
    ts = rand(Nonultrametric(10), 10)
    @test_nowarn show(a, ts)
    @test_nowarn show(a, parsenewick("((,),(,,));", NamedPolytomousTree))
end
end
