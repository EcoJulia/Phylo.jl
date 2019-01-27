module TestIterators

using Phylo
using Compat.Test
using Compat

@testset "$TreeType" for TreeType in [Nonultrametric, Ultrametric]
    ntip = 10
    nur = rand(TreeType(ntip))

    nni = nodenameiter(nur)
    @test Compat.IteratorEltype(nni) == Base.HasEltype()
    @test eltype(nni) == nodenametype(typeof(nur))
    @test Compat.IteratorSize(nni) == Base.HasLength()
    @test length(nni) == 2 * ntip - 1

    ni = nodeiter(nur)
    @test Compat.IteratorEltype(ni) == Base.HasEltype()
    @test eltype(ni) == nodetype(typeof(nur))
    @test Compat.IteratorSize(ni) == Base.HasLength()
    @test length(ni) == 2 * ntip - 1
    @test Set(ni) == Set(map(n -> getnode(nur, n), nni))

    bni = branchnameiter(nur)
    @test Compat.IteratorEltype(bni) == Base.HasEltype()
    @test eltype(bni) == branchnametype(typeof(nur))
    @test Compat.IteratorSize(bni) == Base.HasLength()
    @test length(bni) == ntip * 2 - 2

    bi = branchiter(nur)
    @test Compat.IteratorEltype(bi) == Base.HasEltype()
    @test eltype(bi) == branchtype(typeof(nur))
    @test Compat.IteratorSize(bi) == Base.HasLength()
    @test length(bi) == 2 * ntip - 2
    @test Set(bi) == Set(map(b -> getbranch(nur, b), bni))
end

@testset "$TreeType with filter" for TreeType in [Nonultrametric, Ultrametric]
    ntip = 10
    nur = rand(TreeType(ntip))

    nni = nodenamefilter(n -> isleaf(nur, n), nur)
    @test Compat.IteratorEltype(nni) == Base.HasEltype()
    @test eltype(nni) == nodenametype(typeof(nur))
    @test Compat.IteratorSize(nni) == Base.HasLength()
    @test length(nni) == ntip
    leaves = collect(nni)
    @test Set(leaves) == Set(getleafnames(nur))

    ni = nodefilter(n -> isleaf(nur, n), nur)
    @test Compat.IteratorEltype(ni) == Base.HasEltype()
    @test eltype(ni) == nodetype(typeof(nur))
    @test Compat.IteratorSize(ni) == Base.HasLength()
    @test length(ni) == ntip
    @test all(n -> isleaf(nur, n), collect(ni))

    bni = branchnamefilter(nur) do branch
        return isleaf(nur, dst(nur, branch))
    end
    @test Compat.IteratorEltype(bni) == Base.HasEltype()
    @test eltype(bni) == branchnametype(typeof(nur))
    @test Compat.IteratorSize(bni) == Base.HasLength()
    @test length(bni) == ntip

    bi = branchfilter(nur) do branch
        return isleaf(nur, dst(nur, branch))
    end
    @test Compat.IteratorEltype(bi) == Base.HasEltype()
    @test eltype(bi) == branchtype(typeof(nur))
    @test Compat.IteratorSize(bi) == Base.HasLength()
    @test length(bi) == ntip
    bileaves = map(b -> dst(bi.tree, b), bi)
    @test Set(bileaves) == Set(getnode.(nur, leaves))
end

end
