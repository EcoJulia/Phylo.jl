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
@testset "Nonultrametric{$TreeType}" for TreeType in [
    NamedTree,
    NamedPolytomousTree,
    Phylo.LTD{OneRoot, Float64},
    RootedTree
]
    nt = rand(Nonultrametric{TreeType}(ntips))
    @test_nowarn show(a, nt)
    @test_nowarn show(a, first(nodeiter(nt)))
    @test_nowarn show(a, first(branchiter(nt)))
    ts = rand(Nonultrametric{TreeType}(10), 10)
    @test_nowarn show(a, ts)
    @testset "$newick" for newick in [
        "((A,B),(,,));",
        "((,),(A[&a=\"aaa\", b=1],,));",
        "((,):3.2,(,,):[&a=\"aaa\", b=1]23.2);"
    ]
        ps = parsenewick(newick, TreeType)
        @test_nowarn show(a, ps)
        @test_nowarn show(a, [ps])
        @test_nowarn show(a, first(nodeiter(ps)))
        @test_nowarn show(a, first(branchiter(ps)))
        @test_nowarn show(a, (tree = ps, node = first(getnodes(ps))))
        @test_nowarn show(a,
                          (tree = ps,
                           node = getparent(ps, first(getleaves(ps)))))
        @test_nowarn show(a, (tree = ps, branch = first(getbranches(ps))))
        b = IOBuffer()
        @test_nowarn show(IOContext(b, :compact => false), ps)
        @test Phylo.outputtree(ps, Phylo.StandardOutput()) == String(take!(b))
        @test_nowarn show(IOContext(b, :compact => true), ps)
        @test Phylo.outputtree(ps, Phylo.CompactOutput()) == String(take!(b))
        @test Phylo.outputnode(ps, first(getnodes(ps)),
                               Phylo.CompactOutput()) isa String
        @test Phylo.outputnode(ps, first(getnodes(ps)),
                               Phylo.StandardOutput()) isa String
        @test Phylo.outputbranch(ps, first(getbranches(ps)),
                                 Phylo.CompactOutput()) isa String
        @test Phylo.outputbranch(ps, first(getbranches(ps)),
                                 Phylo.StandardOutput()) isa String
        if TreeType == RootedTree
            @test_nowarn show(a, getnodes(ps))
            @test_nowarn show(a, getbranches(ps))
        end
    end
end
end
