module TestTrim

using Phylo
using DataFrames
using IterableTables: getiterator
using Compat.Test

species = ["Dog", "Cat", "Human"]
ntips = 10
df = DataFrame(species = species, count=[10, 20, 3])
@testset "getinternalnodes()" begin
    # Create a 10 tip tree
    test_tree = rand(Ultrametric(10))
    ints = Phylo.getinternalnodes(test_tree)
    tips = getleafnames(test_tree)
    root = collect(nodenamefilter(isroot, test_tree))
    @test length(ints ∩ tips) == 0
    @test length(ints ∩ root) == 0
end

@testset "keep/droptips!()" begin
    # Create a 10 tip tree
    test_tree = rand(Ultrametric(10))
    keep_tips = ["tip 5", "tip 6", "tip 7", "tip 8", "tip 9", "tip 10"]
    drop_tips = ["tip 1", "tip 2", "tip 3", "tip 4"]
    # Drop all but first three tips
    @test Set(drop_tips) == Set(keeptips!(test_tree, keep_tips))
    tips = getleafnames(test_tree)
    test_tree2 = rand(Ultrametric(10))
    @test Set(drop_tips) == Set(droptips!(test_tree2, drop_tips))
    tips2 = getleafnames(test_tree2)
    @test Set(tips) == Set(keep_tips)
    @test Set(tips) == Set(tips2)
    @test validate(test_tree)
    @test validate(test_tree2)

    tdf = rand(Ultrametric(df))
    @test ["Dog"] == droptips!(tdf, ["Dog"])
    @test length(getiterator(getleafinfo(tdf))) == 2
end

end
