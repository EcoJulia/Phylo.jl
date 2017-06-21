module TestInterface

using Phylo
using Base.Test
using Compat

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
                itr = Compat.Iterators.filter(name ->
                                              hasoutboundspace(tree, name) &&
                                              name != node, allnodes)
                addbranch!(tree, first(itr), node)
            end
        @test Set(branches) == Set(BranchNameIterator(tree))
        @test validate(tree)
        @test_throws ErrorException addbranch!(tree, allnodes[1], allnodes[2])
        who = getparent(tree, species[1])
        b = getinbound(tree, species[1])
        node1 = getnode(tree, species[1])
        @test !isunattached(node1) &&
            hasoutboundspace(node1) && !hasinboundspace(node1) &&
            outdegree(node1) == 0 && indegree(node1) == 1
        @test hasbranch(tree, b)
        @test gettarget(getbranch(tree, b)) == species[1]
        @test getsource(getbranch(tree, b)) == who
        @test b == deletebranch!(tree, b)
        @test !hasbranch(tree, b)
        @test isunattached(tree, species[1]) &&
            hasoutboundspace(tree, species[1]) &&
            hasinboundspace(tree, species[1]) &&
            outdegree(tree, species[1]) == 0 &&
            indegree(tree, species[1]) == 0
        @test_throws ErrorException getbranch(tree, b)
        branches = collect(Compat.Iterators.filter(name -> name != b, branches))
        @test Set(branches) == Set(BranchNameIterator(tree))
        b3 = getinbound(tree, species[2])
        src = getsource(tree, b3)
        tgt = gettarget(tree, b3)
        @test b3 == changesource!(tree, b3, who)
        @test b3 == changetarget!(tree, b3, species[1])
        b2 = addbranch!(tree, src, tgt)
        push!(branches, b2)
        @test Set(branches) == Set(BranchNameIterator(tree))
        @test species[1] == deletenode!(tree, species[1])
        @test !hasnode(tree, species[1])
        @test species[1] == branch!(tree, who, target=species[1])
        @test hasnode(tree, species[1])
        @test validate(tree)
        @test all(map(node -> isleaf(tree, node), species))
        @test all(map(node -> !isroot(tree, node) &&
                      !isunattached(tree, node) &&
                      !isinternal(tree, node), species))
        @test Set(map(pair -> getsource(pair[2]), collect(getbranches(tree)))) ∪
            Set(map(pair -> gettarget(pair[2]), collect(getbranches(tree)))) ==
            Set(allnodes)
    end
end

end
