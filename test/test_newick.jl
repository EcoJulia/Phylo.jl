module TestNewick

using Phylo
using Compat.Test

@testset "newick and nexus parsing" begin
    @testset "For $TreeType" for TreeType in [NamedBinaryTree]
        @test length(nodeiter(parsenewick("((,),(,));"))) == 7
        @test ["where", "when it's good", "Not mine", "MyLeaf"] ⊆
            nodenameiter(parsenewick("""((MyLeaf,"when it's good"),
                                             ('Not mine',where));""",
                                     TreeType))
        tree =
            parsenewick("((MyLeaf[&Real=23,'Not real'={5,4}]:4.0,)Parent,(,));",
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
        if VERSION < v"0.7.0-"
            @test_warn "Tree ended but not finished" parsenewick("((,),(,));a")
        end
        @test_throws Exception parsenewick("((,),(,)));", TreeType)
        @test_throws Exception parsenewick("((,),(,))", TreeType)
        @test_throws Exception parsenewick("((,),(,);", TreeType)
        @test_throws Exception parsenewick("((MyLeaf:-4.0,)Parent,(,));",
                                           TreeType)
        tree = open(f -> parsenewick(f, TreeType), Phylo.path("H1N1.newick"))
        @test nleaves(tree) == 507
        @test ntrees(tree) == 1
        treex = open(f -> parsenexus(f, TreeType), Phylo.path("H1N1.trees"))
        @test nleaves(tree) == nleaves(treex) ==
            nleaves(treex["TREE1"]) == nleaves(treex["TREE2"])
        @test ntrees(treex) == 2
    end
    
    @testset "For $TreeType" for TreeType in
        [NamedTree, RootedTree, ManyRootTree]
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
        @test nleaves(tree) == 507
        @test ntrees(tree) == 1
        treex = open(f -> parsenexus(f, TreeType), Phylo.path("H1N1.trees"))
        @test nleaves(tree) == nleaves(treex) ==
            nleaves(treex["TREE1"]) == nleaves(treex["TREE2"])
        @test ntrees(treex) == 2
    end
end

end

