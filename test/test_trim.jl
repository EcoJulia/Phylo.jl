module TestTrim

using Phylo
using Compat.Test

@testset "getinternalnodes()" begin
    # Create a 10 tip tree
    test_tree = rand(Ultrametric(10))
    ints = getinternalnodes(test_tree)
    tips = getleafnames(test_tree)
    root = collect(nodenamefilter(isroot, test_tree))
    @test length(ints ∩ tips) == 0
    @test length(ints ∩ root) == 0 
end

@testset "drop_tip!()" begin
    # Create a 10 tip tree
    test_tree = rand(Ultrametric(10))
    keep_tips = ["tip 1", "tip 2", "tip 3"]
    # Drop all but first three tips
    drop_tip!(test_tree, keep_tips)
    tips = getleafnames(test_tree)
    @test Set(tips) == Set(keep_tips)
    @test validate(test_tree)
end
