module TestLengths

using Phylo
using Compat.Test

@testset "routes" begin
    species = ["Dog", "Cat", "Human", "Potato"]
    tree = RootedTree(species)
    nr = createnode!(tree)
    n2 = createnode!(tree)
    createbranch!(tree, nr, n2)
    n3 = createnode!(tree)
    createbranch!(tree, n2, n3)
    nh = nodehistory(tree, n2)
    @test nr ∈ nh
    @test Set(nh) == Set(push!(getancestors(tree, n2), n2))
    @test nr ≠ getparent(tree, n3)
    @test n2 == getparent(tree, n3)
    @test nr ∉ getdescendants(tree, nr)
    @test n3 ∈ getdescendants(tree, nr)
    @test n3 ∉ getchildren(tree, nr)
    @test Set([nr, n2, n3]) == Set(nodehistory(tree, n3))
    br = getbranchname(tree, createbranch!(tree, n3, "Dog"))
    @test br ∈ branchhistory(tree, "Dog")
    br2 = getbranchname(tree, createbranch!(tree, n3, "Cat"))
    @test [br2, br] == branchroute(tree, "Cat", "Dog")
    br3 = createbranch!(tree, n2, "Human")
    @test getbranchname(tree, createbranch!(tree, nr, "Potato", name = 551)) == 551
    @test 551 ∈ branchroute(tree, "Dog", "Potato")
    deletebranch!(tree, 551)
    createbranch!(tree, nr, "Potato")
    @test Set(species) ⊆ Set(getnodename.(tree, getdescendants(tree, nr)))
    @test validate!(tree)
end

end
