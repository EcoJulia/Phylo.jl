module TestMetrics

using Phylo
using Test

@testset "Metrics" begin
    tree = open(parsenewick, Phylo.path("H1N1.newick"))
    species = getleaves(tree)[[2, 5, 8, 12, 22]];
    mrca = common_ancestor(tree, species)
    @test getnodename(tree, mrca) == "Node 65"
end

end
