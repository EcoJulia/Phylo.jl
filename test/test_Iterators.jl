module TestIterators

using Phylo
using Base.Test

@testset "$TreeType" for TreeType in [Nonultrametric, Ultrametric]
    ntip = 10
    nur = rand(TreeType(ntip))
    
    nni = NodeNameIterator(nur)
    @test Base.iteratoreltype(nni) == Base.HasEltype()
    @test eltype(nni) == nodenametype(nur)
    @test Base.iteratorsize(nni) == Base.HasLength()
    @test length(nni) == 2 * ntip - 1
    
    ni = NodeIterator(nur)
    @test Base.iteratoreltype(ni) == Base.HasEltype()
    @test eltype(ni) == nodetype(nur)
    @test Base.iteratorsize(ni) == Base.HasLength()
    @test length(ni) == 2 * ntip - 1
    @test Set(ni) == Set(map(n -> getnode(nur, n), nni))

    bni = BranchNameIterator(nur)
    @test Base.iteratoreltype(bni) == Base.HasEltype()
    @test eltype(bni) == branchnametype(nur)
    @test Base.iteratorsize(bni) == Base.HasLength()
    @test length(bni) == ntip * 2 - 2
    
    bi = BranchIterator(nur)
    @test Base.iteratoreltype(bi) == Base.HasEltype()
    @test eltype(bi) == branchtype(nur)
    @test Base.iteratorsize(bi) == Base.HasLength()
    @test length(bi) == 2 * ntip - 2
    @test Set(bi) == Set(map(b -> getbranch(nur, b), bni))
end

@testset "$TreeType with filter" for TreeType in [Nonultrametric, Ultrametric]
    ntip = 10
    nur = rand(TreeType(ntip))
    
    nni = NodeNameIterator(nur, isleaf)
    @test Base.iteratoreltype(nni) == Base.HasEltype()
    @test eltype(nni) == nodenametype(nur)
    @test Base.iteratorsize(nni) == Base.HasLength()
    @test length(nni) == ntip
    leaves = collect(nni)
    @test leaves == getleafnames(nur)
    
    ni = NodeIterator(nur, isleaf)
    @test Base.iteratoreltype(ni) == Base.HasEltype()
    @test eltype(ni) == nodetype(nur)
    @test Base.iteratorsize(ni) == Base.HasLength()
    @test length(ni) == ntip
    @test all(isleaf, collect(ni))

    bni = BranchNameIterator(nur, branch -> isleaf(bni.tree, dst(branch)))
    @test Base.iteratoreltype(bni) == Base.HasEltype()
    @test eltype(bni) == branchnametype(nur)
    @test Base.iteratorsize(bni) == Base.HasLength()
    @test length(bni) == ntip
    
    bi = BranchIterator(nur, branch -> isleaf(bi.tree, dst(branch)))
    @test Base.iteratoreltype(bi) == Base.HasEltype()
    @test eltype(bi) == branchtype(nur)
    @test Base.iteratorsize(bi) == Base.HasLength()
    @test length(bi) == ntip
    bileaves = map(dst, bi)
    @test Set(bileaves) == Set(leaves)
end

end
