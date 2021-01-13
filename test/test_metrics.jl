module TestMetrics

using Phylo
using Test

@testset "Metrics" begin
    trees = open(parsenexus, Phylo.path("H1N1.trees"))
    species = ["H1N1_A_PUERTORICO_8_1934", "H1N1_A_BRAZIL_11_1978"];
    mrca1 = common_ancestor(trees["TREE1"], species)

    # MRCA should be the root
    @test mrca1 == getnodename(tree["TREE1"], getroot(tree["TREE1"]))
    nodes = getnodenames(tree["TREE1"])
    desc1 = getdescendants(tree["TREE1"], mrca1)

    # MRCA + descendants is whole tree
    @test all(sort(nodes) .== sort(desc1 ∪ [mrca1]))
    leaves = getleafnames(tree["TREE1"])
    leaves1 = sort(leaves ∩ desc1)

    # All leaves are in descendants
    @test length(leaves1) == length(leaves)

    # All leaves are the same in both trees
    mrca2 = common_ancestor(trees["TREE2"], reverse(species))
    leaves2 = sort(leaves ∩ getdescendants(tree["TREE2"], mrca2))
    @test all(leaves1 .== leaves2)

    # MRCA for Sendai strains has 8 descendants, including Fukuoka
    desendais =
        getdescendants(tree["TREE1"],
                       common_ancestor(tree["TREE1"],
                                       filter(x -> contains(x, "SENDAI"),
                                       leaves)))
    @test length(desendais) == 8
    @test length(filter(x -> contains(x, "FUKUOKA"), desendais)) == 1
end

end
