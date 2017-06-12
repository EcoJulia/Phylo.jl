module TestRand

using Phylo
using Base.Test

@testset "rand()" begin
    # Create a 10 tip tree
    nu = Nonultrametric(10)
    @test validate(rand(nu))
    @test Set(getleafnames(rand(nu))) == Set(getleafnames(rand(nu)))
    # Create a tree with named tips
    species = ["Dog", "Cat", "Human"]
    t = rand(Nonultrametric(species))
    @test validate(t)
    @test Set(getleafnames(t)) == Set(species)
end

end
