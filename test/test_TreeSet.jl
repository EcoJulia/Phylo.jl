module TestTreeSet

using Phylo
using Compat.Test

@testset "TestSet" begin
    @test length(rand(Ultrametric(10), 10)) ==
        length(rand(Nonultrametric(10), 10))
    ts = rand(Ultrametric(10), 1:10)
    for name in treenameiter(ts)
        @test length(getleafnames(ts[name])) == 10
    end
    for tree in treeiter(ts)
        @test length(getleafnames(tree)) == 10
    end
end
end
