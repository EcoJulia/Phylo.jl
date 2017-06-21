module TestInterface

using Phylo
using Base.Test

@testset "Build and tear down trees" begin
    @testset "For $TreeType" for TreeType in
        [NamedTree,
         BinaryTree{LeafInfo, Nullable{Float64}},
         BinaryTree{LeafInfo, Vector{String}}]

        species = ["Dog", "Cat", "Human", "Potato", "Apple"]
        tree = TreeType(species)
        othernodes = ["Root", "Internal 1", "Internal 2"]
        @test othernodes == addnodes!(tree, othernodes)
        extra = addnodes!(tree, 1)
        @test isa(extra[1], String)
        append!(othernodes, extra)
        also = addnode!(tree)
        @test also ∉ othernodes
        @test_throws ErrorException addnode!(tree, also)
        @test_throws ErrorException addnodes!(tree, [also])
        @test_throws ErrorException deletenode!(tree, "really not")
        @test also == deletenode!(tree, also)
        also = addnode!(tree)
        push!(othernodes, also)
        @test Set(othernodes) ⊆ Set(NodeNameIterator(tree))
        innodes = append!(copy(species), othernodes)
        allnodes = reverse(innodes)
        pop!(innodes)
        branches =
            map(innodes) do node
                itr = Iterators.filter(name -> hasoutboundspace(tree, name) &&
                                       name != node, allnodes)
                addbranch!(tree, first(itr), node)
            end
        @test Set(branches) == Set(BranchNameIterator(tree))
        @test validate(tree)
    end
end

end
