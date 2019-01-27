module TestInterface

using Phylo
using DataFrames
using Compat.Test

@testset "Build and tear down trees" begin
    @testset "For $TreeType" for TreeType in
        [NamedTree, NamedBinaryTree,
         BinaryTree{ManyRoots, DataFrame, Vector{Float64}},
         RootedTree, ManyRootTree]

        species = ["Dog", "Cat", "Human", "Potato", "Apple"]
        tree = TreeType(species)
        othernodes = ["Some 1", "Some 2"]
        nodes = createnodes!(tree, othernodes)
        @test [getnodename(tree, node) for node in nodes] == othernodes
        extra = [getnodename(tree, node) for node in createnodes!(tree, 1)]
        @test isa(extra[1], String)
        append!(othernodes, extra)
        also = getnodename(tree, createnode!(tree))
        @test also ∉ othernodes
        @test_throws ErrorException createnode!(tree, also)
        @test_throws ErrorException createnodes!(tree, [also])
        @test_throws ErrorException deletenode!(tree, "really not")
        @test deletenode!(tree, also)
        also = getnodename(tree, createnode!(tree))
        push!(othernodes, also)
        @test Set(othernodes) ⊆ Set(getnodenames(tree))
        innodes = append!(copy(species), othernodes)
        allnodes = reverse(innodes)
        pop!(innodes)
        @test_throws AssertionError getroot(tree)
        branches = map(innodes) do node
            itr = Iterators.filter(name -> hasoutboundspace(tree, name) &&
                                   name != node, allnodes)
            getbranch(tree, createbranch!(tree, first(itr), node))
        end
        @test_nowarn getroot(tree)
        branchnames = [getbranchname(tree, branch) for branch in branches]
        @test Set(branchnames) == Set(getbranchnames(tree))
        @test validate!(tree)
        @test_throws ErrorException createbranch!(tree, allnodes[1], allnodes[2])
        who = getparent(tree, species[1])
        b = getbranch(tree, getinbound(tree, species[1]))
        bn = getbranchname(tree, b)
        node1 = getnode(tree, species[1])
        @test !isunattached(tree, node1) &&
            hasoutboundspace(tree, node1) && !hasinboundspace(tree, node1) &&
            outdegree(tree, node1) == 0 && indegree(tree, node1) == 1
        @test hasbranch(tree, b)
        @test hasbranch(tree, bn)
        @test getnodename(tree, dst(tree, getbranch(tree, b))) == species[1]
        @test getnodename(tree, src(tree, getbranch(tree, b))) == who
        @test deletebranch!(tree, b)
        @test !hasbranch(tree, b)
        @test isunattached(tree, species[1]) &&
            hasoutboundspace(tree, species[1]) &&
            hasinboundspace(tree, species[1]) &&
            outdegree(tree, species[1]) == 0 &&
            indegree(tree, species[1]) == 0
        @test_throws ErrorException getbranch(tree, b)
        branches = [branch for branch in branches if branch != b]
        @test Set(branches) == Set(getbranches(tree))
        b3 = getinbound(tree, species[2])
        source = src(tree, b3)
        destination = dst(tree, b3)
        @test deletebranch!(tree, b3)
        branches = [branch for branch in branches if branch != b3]
        b3 = createbranch!(tree, who, species[1])
        @test who == getnodename(tree, src(tree, b3))
        spparent = getparent(tree, getnode(tree, species[1]))
        @test getnode(tree, spparent) == getnode(tree, src(tree, b3))
        @test species[1] == getnodename(tree, dst(tree, b3))
        b2 = getbranchname(tree, createbranch!(tree, source, destination))
        branchnames = [getbranchname(tree, branch) for branch in branches]
        push!(branchnames, b2, getbranchname(tree, b3))
        @test Set(branchnames) == Set(getbranchnames(tree))
        @test deletenode!(tree, species[1])
        @test !hasnode(tree, species[1])
        @test species[1] == getnodename(tree, createnode!(tree, species[1]))
        @test_nowarn createbranch!(tree, who, species[1])
        @test hasnode(tree, species[1])
        @test validate!(tree)
        @test all(isleaf(tree, node) for node in species)
        @test all((!isroot(tree, node) & !isunattached(tree, node) &
                   !isinternal(tree, node)) for node in species)
        @test Set(getnodename(tree, src(tree, branch))
                  for branch in getbranches(tree)) ∪
            Set(getnodename(tree, dst(tree, branch))
                for branch in getbranches(tree)) == Set(getnodenames(tree))
        createnode!(tree)
        @test (roottype(TreeType) == OneRoot) ⊻ validate!(tree)
    end
end

end
