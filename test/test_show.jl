module TestShow

using Phylo
using DataFrames

using Test

@info """
These tests only check that the show() commands
do not give warnings or errors, not that they are correct!
"""

ntips = 10
a = IOBuffer()
@testset "Nonultrametric{$TreeType}" for TreeType in
    [NamedTree, NamedPolytomousTree, RootedTree]

    nt = rand(Nonultrametric{TreeType}(ntips))
    @test_nowarn show(a, nt)
    @test_nowarn show(a, first(nodeiter(nt)))
    @test_nowarn show(a, first(branchiter(nt)))
    ts = rand(Nonultrametric{TreeType}(10), 10)
    @test_nowarn show(a, ts)
    ps = parsenewick("((,),(,,));", TreeType)
    @test_nowarn show(a, ps)
    @test_nowarn show(a, first(nodeiter(ps)))
    @test_nowarn show(a, first(branchiter(ps)))
end
end
