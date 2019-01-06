module TestInterface

using Phylo
using DataFrames
using Compat.Test

@testset "Build and tear down trees" begin
    @testset "For $TreeType" for TreeType in
        [NamedTree, NamedBinaryTree,
         BinaryTree{ManyRoots, DataFrame, Vector{Float64}}]

        species = ["Dog", "Cat", "Human", "Potato", "Apple"]
        tree = TreeType(species)
        othernodes = ["Some 1", "Some 2"]
        @test othernodes == createnodes!(tree, othernodes)
        extra = createnodes!(tree, 1)
        @test isa(extra[1], String)
        append!(othernodes, extra)
        also = createnode!(tree)
        @test also ∉ othernodes
        @test_throws ErrorException createnode!(tree, also)
        @test_throws ErrorException createnodes!(tree, [also])
        @test_throws ErrorException deletenode!(tree, "really not")
        @test also == deletenode!(tree, also)
        also = createnode!(tree)
        push!(othernodes, also)
        @test Set(othernodes) ⊆ Set(getnodenames(tree))
        innodes = append!(copy(species), othernodes)
        allnodes = reverse(innodes)
        pop!(innodes)
        branches =
            map(innodes) do node
                itr = Iterators.filter(name ->
                                       hasoutboundspace(tree, name) &&
                                       name != node, allnodes)
                createbranch!(tree, first(itr), node)
            end
        @test Set(branches) == Set(getbranchnames(tree))
        @test validate(tree)
        @test_throws ErrorException createbranch!(tree, allnodes[1], allnodes[2])
        who = getparent(tree, species[1])
        b = getinbound(tree, species[1])
        node1 = getnode(tree, species[1])
        @test !isunattached(tree, node1) &&
            hasoutboundspace(tree, node1) && !hasinboundspace(tree, node1) &&
            outdegree(tree, node1) == 0 && indegree(tree, node1) == 1
        @test hasbranch(tree, b)
        @test dst(getbranch(tree, b)) == species[1]
        @test src(getbranch(tree, b)) == who
        @test b == deletebranch!(tree, b)
        @test !hasbranch(tree, b)
        @test isunattached(tree, species[1]) &&
            hasoutboundspace(tree, species[1]) &&
            hasinboundspace(tree, species[1]) &&
            outdegree(tree, species[1]) == 0 &&
            indegree(tree, species[1]) == 0
        @test_throws ErrorException getbranch(tree, b)
        branches = [name for name in branches if name != b]
        @test Set(branches) == Set(getbranches(tree))
        b3 = getinbound(tree, species[2])
        source = src(tree, b3)
        destination = dst(tree, b3)
        @test b3 == deletebranch!(tree, b3)
        b3 = createbranch!(tree, who, species[1])
        @test who == src(tree, b3)
        @test species[1] == dst(tree, b3)
        b2 = createbranch!(tree, source, destination)
        push!(branches, b2)
        @test Set(branches) == Set(getbranchnames(tree))
        @test species[1] == deletenode!(tree, species[1])
        @test !hasnode(tree, species[1])
        @test species[1] == createnode!(tree, species[1])
        @test createbranch!(tree, who, species[1]) isa Int
        @test hasnode(tree, species[1])
        @test validate(tree)
        @test all(isleaf(tree, node) for node in species)
        @test all((!isroot(tree, node) & !isunattached(tree, node) &
                   !isinternal(tree, node)) for node in species)
        @test Set(src(tree, branch) for branch in getbranches(tree)) ∪
            Set(dst(tree, branch) for branch in getbranches(tree)) ==
            Set(allnodes)
    end
end

end
