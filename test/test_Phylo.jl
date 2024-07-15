# SPDX-License-Identifier: BSD-2-Clause

module TestPhylo
using Test

using Phylo

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

end
