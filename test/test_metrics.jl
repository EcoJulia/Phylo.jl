module TestMetrics

using Phylo
using Test
using Compat

@testset "Metrics" begin
    trees = open(parsenexus, Phylo.path("H1N1.trees"))
    tree1, tree2 = trees["TREE1"], trees["TREE2"]
    species = ["H1N1_A_PUERTORICO_8_1934", "H1N1_A_BRAZIL_11_1978"];
    mrca1 = mrca(tree1, species)

    # MRCA should be the root
    @test mrca1 == getnodename(tree1, getroot(tree1))
    nodes = getnodenames(tree1)
    desc1 = getdescendants(tree1, mrca1)

    # MRCA + descendants is whole tree
    @test all(sort(nodes) .== sort(desc1 ∪ [mrca1]))
    leaves = getleafnames(tree1)
    leaves1 = sort(leaves ∩ desc1)

    # All leaves are in descendants
    @test length(leaves1) == length(leaves)

    # All leaves are the same in both trees
    mrca2 = mrca(tree2, reverse(species))
    leaves2 = sort(leaves ∩ getdescendants(tree2, mrca2))
    @test all(leaves1 .== leaves2)

    # MRCA for Sendai strains has 8 descendants, including Fukuoka
    desendais =
        getdescendants(tree1,
                       mrca(tree1,
                                       filter(x -> Compat.contains(x, "SENDAI"),
                                       leaves)))
    @test length(desendais) == 8
    @test length(filter(x -> Compat.contains(x, "FUKUOKA"), desendais)) == 1
end

end
