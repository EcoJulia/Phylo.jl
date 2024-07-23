# SPDX-License-Identifier: BSD-2-Clause

module TestPhylo
using Test
using Phylo
using Random

@testset "Deprecations" begin
    t = NamedTree()
    @test_deprecated name = addnode!(t)
    @test_deprecated addnode!(t, "Second")
    @test_deprecated addbranch!(t, name, "Second")
    @test_deprecated branch!(t, name)
    @test nleaves(t) == 2
    @test nroots(t) == 1
    @test nbranches(t) == 2
    @test_deprecated treeiter(t)
    @test_deprecated treenameiter(t)
end

@testset "Sort $TreeType" for TreeType in [RootedTree, Phylo.LTD{OneRoot, Float64}]
    tree = rand(Nonultrametric{TreeType}(100))
    @test validate!(deepcopy(tree))
    t2 = sort(tree)
    leaves = getleaves(t2, inorder)
    @test maximum(getheight.(Ref(t2), leaves)) == getheight(t2, leaves[end])
    t3 = sort(tree, uselength = false, rev = true)
    leaves = getleaves(t3, inorder)
    for b in getbranches(t3)
        b.length = one(b.length)
    end
    @test maximum(getheight.(Ref(t3), leaves)) == getheight(t3, first(leaves))
end

end
