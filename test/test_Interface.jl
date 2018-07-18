module TestInterface

using Phylo
using DataFrames
using Compat.Test

@testset "Build and tear down trees" begin
    @testset "For $TreeType" for TreeType in
        [NamedTree,
         BinaryTree{DataFrame, Vector{Float64}},
         BinaryTree{DataFrame, Vector{String}}]

        species = ["Dog", "Cat", "Human", "Potato", "Apple"]
        tree = TreeType(species)
        othernodes = ["Some 1", "Some 2"]
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
        @test Set(othernodes) ⊆ Set(nodenameiter(tree))
        innodes = append!(copy(species), othernodes)
        allnodes = reverse(innodes)
        pop!(innodes)
        branches =
            map(innodes) do node
                itr = Iterators.filter(name ->
                                       hasoutboundspace(tree, name) &&
                                       name != node, allnodes)
                addbranch!(tree, first(itr), node)
            end
        @test Set(branches) == Set(branchnameiter(tree))
        @test validate(tree)
        @test_throws ErrorException addbranch!(tree, allnodes[1], allnodes[2])
        who = getparent(tree, species[1])
        b = getinbound(tree, species[1])
        node1 = getnode(tree, species[1])
        @test !isunattached(node1) &&
            hasoutboundspace(node1) && !hasinboundspace(node1) &&
            outdegree(node1) == 0 && indegree(node1) == 1
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
        branches = collect(Iterators.filter(name -> name != b, branches))
        @test Set(branches) == Set(branchnameiter(tree))
        b3 = getinbound(tree, species[2])
        source = src(tree, b3)
        destination = dst(tree, b3)
        @test b3 == changesrc!(tree, b3, who)
        @test b3 == changedst!(tree, b3, species[1])
        b2 = addbranch!(tree, source, destination)
        push!(branches, b2)
        @test Set(branches) == Set(branchnameiter(tree))
        @test species[1] == deletenode!(tree, species[1])
        @test !hasnode(tree, species[1])
        @test species[1] == branch!(tree, who, destination=species[1])
        @test hasnode(tree, species[1])
        @test validate(tree)
        @test all(map(==, Pair(first(branchiter(tree))),
                      Tuple(tree, first(branchnameiter(tree)))))
        @test all(map(==, Pair(tree, first(branchnameiter(tree))),
                      Tuple(first(branchiter(tree)))))
        @test all(map(node -> isleaf(tree, node), species))
        @test all(map(node -> !isroot(tree, node) &&
                      !isunattached(tree, node) &&
                      !isinternal(tree, node), species))
        @test Set(map(pair -> src(pair[2]), collect(getbranches(tree)))) ∪
            Set(map(pair -> dst(pair[2]), collect(getbranches(tree)))) ==
            Set(allnodes)
    end
end

end
