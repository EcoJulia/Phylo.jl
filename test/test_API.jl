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

@testset "Check errors" begin
    tt = TestTree();
    tn = TestNode();
    tb = TestBranch();
    
    @test_throws ErrorException _getnodes(tt)
    @test_throws ErrorException _getnodenames(tt)
    @test_throws ErrorException _getbranches(tb)
    @test_throws ErrorException _getnodes(tt)
    @test_throws ErrorException _hasinbound(tt, tn)
    @test_throws ErrorException _getinbound(tt, tn)
    @test_throws ErrorException _getoutbounds(tt, tn)
    @test_throws ErrorException _addconnection!(tt, tn, tn)
    @test_throws ErrorException _removeconnection!(tt, tn, tb)
    @test_throws ErrorException _src(tt, tb)
    @test_throws ErrorException _dst(tt, tb)
    @test_throws ErrorException _getnodedata(tt, tn)
    @test_throws ErrorException _getbranchdata(tt, tb)
end

end