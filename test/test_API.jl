module TestAPI

using Phylo
using Phylo.API
using Test

struct TestBranch <: Phylo.AbstractBranch{OneRoot, String}
end

struct TestNode <: Phylo.AbstractNode{OneRoot, String}
end

struct TestTree <: Phylo.AbstractTree{OneTree, OneRoot, String, TestNode, TestBranch}
end

import Phylo.API: _preferbranchobjects
_preferbranchobjects(::Type{<: TestBranch}) = false

@testset "Check errors" begin
    tt = TestTree();
    tn = TestNode();
    tb = TestBranch();
    
    @test_throws MethodError _getnodes(tt, preorder)
    @test_throws ErrorException _getnodenames(tt)
    @test_throws MethodError _nnodes(tt)
    @test_throws ErrorException _getbranches(tt)
    @test_throws ErrorException _getbranchnames(tt)
    @test_throws ErrorException _nbranches(tt)
    @test_throws MethodError _hasnode(tt, tn)
    @test_throws ErrorException _hasnode(tt, "test")
    @test_throws MethodError _hasbranch(tt, tb)
    @test_throws ErrorException _hasbranch(tt, 1)
    @test_throws ErrorException _hasinbound(tt, tn)
    @test_throws ErrorException _getinbound(tt, tn)
    @test_throws ErrorException _getoutbounds(tt, tn)
    @test_throws ErrorException _addconnection!(tt, tn, tn)
    @test_throws ErrorException _removeconnection!(tt, tn, tb)
    @test_throws ErrorException _src(tt, tb)
    @test_throws ErrorException _dst(tt, tb)
    @test_throws ErrorException _getnodedata(tt, tn)
    @test_throws ErrorException _getbranchdata(tt, tb)
    @test !_renamenode!(tt, tn, "New")
end

end