# SPDX-License-Identifier: BSD-2-Clause

module TestMetrics

using Phylo
using AxisArrays: axisvalues
using Test

@testset "Metrics" begin
    trees = open(parsenexus, Phylo.path("H1N1.trees"))
    tree1, tree2 = trees["TREE1"], trees["TREE2"]
    species = ["H1N1_A_PUERTORICO_8_1934", "H1N1_A_BRAZIL_11_1978"]
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
    desendais = getdescendants(tree1,
                               mrca(tree1,
                                    filter(x -> occursin("SENDAI", x),
                                           leaves)))
    @test length(desendais) == 8
    @test length(filter(x -> occursin("FUKUOKA", x), desendais)) == 1
end

# Distances summed along the route between each pair of leaves, which is how
# distances() computes them for trees without the fast path
function routedistances(tree)
    leaves = getleaves(tree)
    return [distance(tree, li, lj) for li in leaves, lj in leaves]
end

@testset "distances and vcv" begin
    trees = open(parsenexus, Phylo.path("H1N1.trees"))
    polytomy = parsenewick("((A:1.0,B:2.0,C:0.5)X:1.0,(D:1.0,E:1.5)Y:2.0)R;")
    singlechild = parsenewick("(((A:1.0,B:2.0):0.4)K:0.6,(C:1.0,D:1.5):2.0)R;")
    for tree in (trees["TREE1"], rand(Nonultrametric(60)),
                 rand(Ultrametric(60)), polytomy, singlechild)
        names = getnodename.(tree, getleaves(tree))
        d = distances(tree)
        V = vcv(tree)
        # Same values, axes and leaf order as summing each route
        @test d ≈ routedistances(tree)
        @test collect(axisvalues(d)[1]) == names
        @test collect(axisvalues(V)[1]) == names
        @test all(iszero, [d[i, i] for i in axes(d, 1)])
        # The diagonal of vcv is the leaf heights, and distances follow from it
        h = [V[i, i] for i in axes(V, 1)]
        leafheights = heightstoroot(tree)
        @test h ≈ [leafheights[name] for name in names]
        @test d ≈ h .+ h' .- 2 .* V
        @test V == V'
        # Off the diagonal, vcv is the height of the most recent common ancestor
        heights = nodeheights(tree)
        for _ in 1:20
            i, j = rand(eachindex(names)), rand(eachindex(names))
            i == j && continue
            @test V[i, j] ≈ heights[mrca(tree, [names[i], names[j]])]
        end
    end

    # Without branch lengths distances() still sums routes, and vcv() refuses
    nolengths = RootedTree(["A", "B", "C"])
    n = createnode!(nolengths)
    createbranch!(nolengths, n, "A")
    createbranch!(nolengths, n, "B")
    r = createnode!(nolengths)
    createbranch!(nolengths, r, n)
    createbranch!(nolengths, r, "C")
    @test isequal(collect(distances(nolengths)), routedistances(nolengths))
    @test_throws ArgumentError vcv(nolengths)
end

end
