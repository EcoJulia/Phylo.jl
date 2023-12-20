module TestNexus

using Phylo
using Test

@testset "Nexus parsing" begin
    @testset "For $TreeType" for TreeType in [NamedBinaryTree, RootedTree]
        tree = open(f -> parsenewick(f, TreeType), Phylo.path("H1N1.newick"))
        @test nleaves(tree) == 507
        @test ntrees(tree) == 1
        treex = open(f -> parsenexus(f, TreeType), Phylo.path("H1N1.trees"))
        @test nleaves(tree) == nleaves(treex) ==
            nleaves(treex["TREE1"]) == nleaves(treex["TREE2"])
        @test ntrees(treex) == 2
        @test Set(gettreenames(treex)) == Set(["TREE1", "TREE2"])
        @test gettreeinfo(treex)["TREE1"]["lnP"] ≈ -gettreeinfo(treex)["TREE2"]["lnP"]
    end
    
    @testset "For $TreeType" for TreeType in
        [NamedTree, Phylo.LTD{OneRoot, Float64}, Phylo.LTD{ManyRoots, Float64}, RootedTree, ManyRootTree]
        tree = open(f -> parsenewick(f, TreeType), Phylo.path("H1N1.newick"))
        tree_p = open(f -> parse(TreeType, f, format = Newick()), Phylo.path("H1N1.newick"))
        @test getnodenames(tree) == getnodenames(tree_p)
        tree_pf = open(parse(TreeType), Phylo.path("H1N1.newick"))
        @test getnodenames(tree) == getnodenames(tree_pf)
        treeset = TreeSet([tree])
        @test gettreename(treeset) == 1
        @test gettree(treeset) ≡ treeset[1]
        treex = open(f -> parsenexus(f, TreeType), Phylo.path("H1N1.trees"))
        @test_throws AssertionError gettreename(treex)
        @test Set(getleafnames(treex)) == Set(getleafnames(treex["TREE1"]))
        treex_p = open(parse(Phylo.treesettype(TreeType), format = Nexus()), Phylo.path("H1N1.trees"))
        @test getnodenames(treex) == getnodenames(treex_p)
        treex_pf = open(f -> parse(Phylo.treesettype(TreeType), f), Phylo.path("H1N1.trees"))
        @test getnodenames(treex) == getnodenames(treex_pf)
        a = IOBuffer()
        @test_nowarn show(a, treex)
        @test ntrees(treex) == 2
        @test nleaves(tree) == nleaves(treex) ==
            nleaves(treex["TREE1"]) == nleaves(treex["TREE2"])
        @test nnodes(treex, "TREE1") == nnodes(treex)["TREE1"]
        @test getnodes(treex, "TREE1") == getnodes(treex)["TREE1"]
        @test Set(getleafnames(treex)) ==
            Set(getleafnames(treex["TREE1"])) == Set(getleafnames(treex["TREE2"]))
        if roottype(TreeType) == OneRoot
            @test_nowarn Phylo.write("t.trees", treex)
            treex2 = open(f -> parsenexus(f, TreeType), "t.trees")
            @test_nowarn show(a, treex)
            @test ntrees(treex2) == 2
            @test nleaves(treex2) == nleaves(treex) ==
                nleaves(treex2["TREE1"]) == nleaves(treex2["TREE2"])
            @test Set(getleafnames(treex)) == Set(getleafnames(treex2)) ==
                Set(getleafnames(treex2["TREE1"])) == Set(getleafnames(treex2["TREE2"]))
            @test gettreeinfo(treex2)["TREE1"]["lnP"] ≈ -gettreeinfo(treex2)["TREE2"]["lnP"]
        end
    end
end

end

