module TestInterface

using Phylo
using DataFrames
using Test

matchbranch = true
@testset "Build and tear down trees" begin
    @testset "For $TreeType" for TreeType in
        [NamedTree, NamedBinaryTree,
         BinaryTree{ManyRoots, DataFrame, Vector{Float64}},
         RootedTree, ManyRootTree]

        @test treetype(TreeType) == OneTree
        species = ["Dog", "Cat", "Human", "Potato", "Apple"]
        tree = TreeType(species)
        othernodes = ["Some 1", "Some 2"]
        nodes = createnodes!(tree, othernodes)
        @test !((roottype(TreeType) ≡ ManyRoots) ⊻ validate!(tree))
        @test [getnodename(tree, node) for node in nodes] == othernodes
        extra = [getnodename(tree, node) for node in createnodes!(tree, 1)]
        @test isa(extra[1], String)
        append!(othernodes, extra)
        also = getnodename(tree, createnode!(tree))
        @test also ∉ othernodes
        @test_throws Exception createnode!(tree, also)
        @test_throws Exception createnodes!(tree, [also])
        @test_throws Exception deletenode!(tree, "really not")
        @test deletenode!(tree, also)
        also = getnodename(tree, createnode!(tree))
        push!(othernodes, also)
        @test Set(othernodes) ⊆ Set(getnodenames(tree))
        innodes = append!(copy(species), othernodes)
        allnodes = copy(innodes)
        pop!(innodes)
        @test_throws Exception getroot(tree)
        branches = map(innodes) do node
            itr = Iterators.filter(name -> hasoutboundspace(tree, name) &&
                                   name != node &&
                                   name ∉ getdescendants(tree, node) &&
                                   name ∉ species, allnodes)
            global matchbranch = !matchbranch
            if matchbranch
                getbranch(tree, createbranch!(tree, getnode(tree, first(itr)), getnode(tree, node)))
            else
                getbranch(tree, createbranch!(tree, getnodename(tree, first(itr)), getnodename(tree, node)))
            end
        end
        @test_nowarn getroot(tree)
        branchnames = [getbranchname(tree, branch) for branch in branches]
        @test Set(branchnames) == Set(getbranchnames(tree))
        @test validate!(tree)
        @test_throws Exception createbranch!(tree, allnodes[1], allnodes[2])
        who = getparent(tree, species[1])
        b = getbranch(tree, getinbound(tree, species[1]))
        bn = getbranchname(tree, b)
        node1 = getnode(tree, species[1])
        @test !isunattached(tree, node1) &&
            hasoutboundspace(tree, node1) && !hasinboundspace(tree, node1) &&
            outdegree(tree, node1) == 0 && indegree(tree, node1) == 1 && degree(tree, node1) == 1
        @test hasbranch(tree, b)
        @test hasbranch(tree, bn)
        @test getnodename(tree, dst(tree, getbranch(tree, b))) == species[1]
        @test src(tree, getbranch(tree, b)) == getnode(tree, who)
        @test deletebranch!(tree, b)
        @test !hasbranch(tree, b)
        @test isunattached(tree, species[1]) &&
            hasoutboundspace(tree, species[1]) &&
            hasinboundspace(tree, species[1]) &&
            outdegree(tree, species[1]) == 0 &&
            indegree(tree, species[1]) == 0 &&
            degree(tree, node1) == 0
        @test_throws Exception getbranch(tree, b)
        branches = [branch for branch in branches if branch != b]
        @test Set(branches) == Set(getbranches(tree))
        b3 = getinbound(tree, species[2])
        source = src(tree, b3)
        destination = dst(tree, b3)
        @test deletebranch!(tree, b3)
        branches = [branch for branch in branches if branch != b3]
        createbranch!(tree, who, species[1])
        deletebranch!(tree, who, species[1])
        @test isunattached(tree, species[1]) &&
            hasoutboundspace(tree, species[1]) &&
            hasinboundspace(tree, species[1]) &&
            outdegree(tree, species[1]) == 0 &&
            indegree(tree, species[1]) == 0 &&
            degree(tree, node1) == 0
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
        NT = typeof(node1)
        @test NT ≡ String || nodenametype(NT) ≡ String
        @test NT ≡ String || branchnametype(NT) ≡ Int
        @test NT ≡ String || roottype(NT) ≡ roottype(TreeType)
        BT = typeof(b)
        @test nodenametype(BT) ≡ String
        @test branchnametype(BT) ≡ Int
        @test roottype(BT) ≡ roottype(TreeType)
        @test all(isleaf(tree, node) for node in species)
        @test all((!isroot(tree, node) & !isunattached(tree, node) &
                   !isinternal(tree, node)) for node in species)
        @test Set(getnodename(tree, src(tree, branch))
                  for branch in getbranches(tree)) ∪
            Set(getnodename(tree, dst(tree, branch))
                for branch in getbranches(tree)) == Set(getnodenames(tree))
        createnode!(tree)
        @test treenametype(TreeType) ≡ Int ? gettreename(tree) == 1 : gettreename(tree) == "Tree"
        @test (roottype(TreeType) ≡ OneRoot) ⊻ validate!(tree)
        @test gettree(tree) ≡ tree
        @test length(getnodes(tree)) == nnodes(tree)
        @test ninternal(tree) == nnodes(tree) - nleaves(tree)
    end
end

end
