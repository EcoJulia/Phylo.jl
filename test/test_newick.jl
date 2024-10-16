# SPDX-License-Identifier: BSD-2-Clause

module TestNewick

using Phylo
using Test

@testset "Newick parsing" begin
    @testset "For $TreeType" for TreeType in [NamedBinaryTree]
        @test length(nodeiter(parsenewick("((,),(,));"))) == 7
        @test ["where", "when it's good", "Not mine", "MyLeaf"] ⊆
              nodenameiter(parsenewick("""((MyLeaf,"when it's good"),
                                               ('Not mine',where));""",
                                       TreeType))
        tree = parsenewick("((MyLeaf[&Real=23,'Not real'={5,4}]:4.0,)Parent,(,));",
                           TreeType)
        branches = branchfilter(tree) do tree, branch
            return getnodename(tree, src(tree, branch)) == "Parent" &&
                   getnodename(tree, dst(tree, branch)) == "MyLeaf"
        end
        @test length(branches) == 1
        @test getlength(tree, first(branches)) ≈ 4.0
        @test getnodedata(tree, "MyLeaf")["Real"] == 23
        @test 5 ∈ getnodedata(tree, "MyLeaf")["Not real"]
        @test 4 ∈ getnodedata(tree, "MyLeaf")["Not real"]
        @test_throws Exception parsenewick("((,),(,)));", TreeType)
        @test_throws Exception parsenewick("((,),(,))", TreeType)
        @test_throws Exception parsenewick("((,),(,);", TreeType)
        @test_throws Exception parsenewick("((MyLeaf:-4.0,)Parent,(,));",
                                           TreeType)
        tree = open(f -> parsenewick(f, TreeType), Phylo.path("H1N1.newick"))
        @test nleaves(tree) == 507
        @test ntrees(tree) == 1
    end

    @testset "For $TreeType" for TreeType in [
        NamedTree,
        Phylo.LTD{OneRoot, Float64},
        Phylo.LTD{ManyRoots, Float64},
        RootedTree,
        ManyRootTree
    ]
        @test nnodes(parsenewick("((,),(,,));", TreeType)) == 8
        @test ["where", "when it's good", "Not mine", "MyLeaf", "next"] ⊆
              nodenameiter(parsenewick("""((MyLeaf,"when it's good",next),
                                               ('Not mine',where));""",
                                       TreeType))
        tree = parsenewick("((MyLeaf:4.0,)Parent,());", TreeType)
        branches = branchfilter(tree) do tree, branch
            return getnodename(tree, src(tree, branch)) == "Parent" &&
                   getnodename(tree, dst(tree, branch)) == "MyLeaf"
        end
        @test length(branches) == 1
        @test getlength(tree, first(branches)) ≈ 4.0
        @test_throws Exception parsenewick("((,),(,)));", TreeType)
        @test_throws Exception parsenewick("((,),(,))", TreeType)
        @test_throws Exception parsenewick("((,),(,);", TreeType)
        tree = open(f -> parsenewick(f, TreeType), Phylo.path("H1N1.newick"))
        a = IOBuffer()
        @test_nowarn show(IOContext(a, :compact => false), tree)
        @test String(take!(a)) == Phylo.outputtree(tree, Phylo.StandardOutput())
        @test_nowarn show(IOContext(a, :compact => true), tree)
        @test String(take!(a)) == Phylo.outputtree(tree, Phylo.CompactOutput())
        @test Phylo.outputnode(tree, first(getroots(tree)),
                               Phylo.CompactOutput()) isa String
        @test Phylo.outputnode(tree, first(getroots(tree)),
                               Phylo.StandardOutput()) isa String
        @test Phylo.outputnode(tree, first(getleaves(tree)),
                               Phylo.CompactOutput()) isa String
        @test Phylo.outputnode(tree, first(getleaves(tree)),
                               Phylo.StandardOutput()) isa String
        @test Phylo.outputbranch(tree, first(getbranches(tree)),
                                 Phylo.CompactOutput()) isa String
        @test Phylo.outputbranch(tree, first(getbranches(tree)),
                                 Phylo.StandardOutput()) isa String
        @test Phylo.outputnode(tree, first(getinternalnodes(tree)),
                               Phylo.CompactOutput()) isa String
        @test Phylo.outputnode(tree, first(getinternalnodes(tree)),
                               Phylo.StandardOutput()) isa String
        if roottype(TreeType) == OneRoot
            @test_nowarn Phylo.write("t1.newick", tree)
            tree2 = open(f -> parsenewick(f, TreeType), "t1.newick")
            @test nleaves(tree) == nleaves(tree2) == 507
            @test ntrees(tree) == ntrees(tree2) == 1
            @test_nowarn Phylo.write("t1.trees", tree, format = Nexus())
            tree2x = open(f -> parsenexus(f, TreeType), "t1.trees")
            @test nleaves(tree) == nleaves(tree2x) == 507
            @test ntrees(tree) == ntrees(tree2x) == 1
            @test Set(getnodenames(tree)) == Set(getnodenames(tree2))
            @test Set(getleafnames(tree)) == Set(getleafnames(tree2)) ==
                  Set(getleafnames(tree2x))
            names = getnodenames(tree2)
            Phylo.write("t2.newick", tree2,
                        format = Newick(Dict(names .=> 1:nnodes(tree2))))
            tree3 = open(f -> parsenewick(f, TreeType), "t2.newick")
            @test Set(getnodenames(tree3)) == Set(string.(1:nnodes(tree3)))
            Phylo.write("t3.newick", tree3,
                        format = Newick(Dict(string.(1:nnodes(tree3)) .=> names)))
            tree4 = open(f -> parsenewick(f, TreeType), "t3.newick")
            @test Set(getnodenames(tree)) == Set(getnodenames(tree4))
            @test Set(getleafnames(tree)) == Set(getleafnames(tree4))
        else
            @test nleaves(tree) == 507
            @test ntrees(tree) == 1
        end
    end
end

end
