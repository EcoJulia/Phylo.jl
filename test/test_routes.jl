# SPDX-License-Identifier: BSD-2-Clause

module TestRoutes

using Phylo
using Test

@testset "Routes" begin
    @testset "For $TreeType" for TreeType in [NamedTree, NamedBinaryTree,
        Phylo.LTD{OneRoot, Float64}, Phylo.LTD{ManyRoots, Float64},
        RootedTree, ManyRootTree]
        species = ["Dog", "Cat", "Human", "Potato"]
        tree = TreeType(species)
        nr = createnode!(tree)
        n2 = createnode!(tree)
        createbranch!(tree, getnodename(tree, nr), getnodename(tree, n2))
        n3 = createnode!(tree)
        createbranch!(tree, getnode(tree, n2), getnode(tree, n3))
        nh = nodehistory(tree, n2)
        @test nr ∈ nh
        @test Set(nh) == Set(push!(getancestors(tree, n2), n2))
        @test nr ≠ getparent(tree, n3)
        @test n2 == getparent(tree, n3)
        @test nr ∉ getdescendants(tree, nr)
        @test n3 ∈ getdescendants(tree, nr)
        @test n3 ∉ getchildren(tree, nr)
        @test Set([nr, n2, n3]) == Set(nodehistory(tree, n3))
        br = getbranch(tree, createbranch!(tree, n3, "Dog"))
        @test br ∈ branchhistory(tree, "Dog")
        br2 = getbranch(tree, createbranch!(tree, n3, "Cat"))
        @test [br2, br] == branchroute(tree, "Cat", "Dog")
        br3 = createbranch!(tree, n2, "Human")
        @test getbranchname(tree,
                            createbranch!(tree, nr, "Potato",
                                          name = 551)) == 551
        @test 551 ∈ getbranchname.(tree, branchroute(tree, "Dog", "Potato"))
        deletebranch!(tree, 551)
        createbranch!(tree, nr, "Potato")
        @test Set(species) ⊆ Set(getnodename.(tree, getdescendants(tree, nr)))
        @test validate!(tree)
        root = getroot(tree)
        @test Set(getdescendants(tree, getchildren(tree, root)[1])) ⊆
              Set(getdescendants(tree, root)) ⊆ Set(traversal(tree))
        @test Set(getdescendants(tree,
                                 getchildren(tree,
                                             getnodename(tree, root))[1])) ⊆
              Set(getdescendants(tree, getnodename(tree, root))) ⊆
              Set(getnodename.(tree, traversal(tree)))
        @test length(collect(traversal(tree))) == nnodes(tree)
        @test length(getdescendants(tree, root)) ==
              length(branchfuture(tree, root))
        nn = first(nodenamefilter(isleaf, tree))
        @test length(branchhistory(tree, nn)) == length(getancestors(tree, nn))
        @test length(branchfuture(tree, root)) ==
              length(nodefuture(tree, root)) - 1
        @test Set(getancestors(tree, nn)) ⊆ Set(nodehistory(tree, nn))
        n = first(nodefilter(isleaf, tree))
        @test length(branchhistory(tree, n)) == length(getancestors(tree, n))
        @test Set(getancestors(tree, n)) ⊆ Set(nodehistory(tree, n))
        @test Set(traversal(tree, anyorder)) ==
              Set(traversal(tree, preorder)) ==
              Set(traversal(tree, inorder)) ==
              Set(traversal(tree, postorder)) ==
              Set(traversal(tree, breadthfirst))
    end
end

end
