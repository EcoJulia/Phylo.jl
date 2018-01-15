module TestLengths

using Phylo
using Compat.Test

@testset "routes" begin
    species = ["Dog", "Cat", "Human", "Potato"]
    tree = NamedTree(species)
    nr = addnode!(tree)
    n2 = branch!(tree, nr)
    n3 = branch!(tree, n2)
    nh = nodehistory(tree, n2)
    @test nr ∈ nh
    @test Set(nh) == Set(push!(getancestors(tree, n2), n2))
    @test nr ≠ getparent(tree, n3)
    @test n2 == getparent(tree, n3)
    @test nr ∉ getdescendants(tree, nr)
    @test n3 ∈ getdescendants(tree, nr)
    @test n3 ∉ getchildren(tree, nr)
    @test Set([nr, n2, n3]) == Set(nodehistory(tree, n3))
    br = addbranch!(tree, n3, "Dog")
    @test br ∈ branchhistory(tree, "Dog")
    br2 = addbranch!(tree, n3, "Cat")
    @test [br2, br] == branchroute(tree, "Cat", "Dog")
    br3 = addbranch!(tree, n2, "Human")
    @test addbranch!(tree, nr, "Potato", branchname = 551) == 551
    @test 551 ∈ branchroute(tree, "Dog", "Potato")
    deletebranch!(tree, 551)
    br4 = addbranch!(tree, nr, "Potato")
    @test Set(species) ⊆ Set(getdescendants(tree, nr))
    @test validate(tree)
end

end
