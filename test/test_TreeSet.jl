module TestTreeSet

using Phylo
using Compat.Test

@testset "TestSet" begin
    @test length(rand(Ultrametric(10), 10)) ==
        length(rand(Nonultrametric(10), 10))
    ts = rand(Ultrametric(10), 1:10)
    @test ntrees(ts) == 10
    @test length(getbranchnames(ts)) == 18
    for name in treenameiter(ts)
        @test name ∈ 1:10
        @test length(getleafnames(ts[name])) == 10
    end
    for tree in treeiter(ts)
        @test length(getnodenames(tree)) == 19
    end
    for tree in ts
        @test length(getbranchnames(tree)) == 18
    end
    @test getleafnames(TreeSet(Dict{String, NamedTree}(), Dict{String, Any}())) == String[]
    @test getbranchnames(TreeSet(Dict{String, NamedTree}(), Dict{String, Any}())) == Int[]
end
end
