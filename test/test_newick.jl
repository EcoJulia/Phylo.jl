module TestNewick
using Base.Test

using Phylo

@testset "A few simple trees" begin
    @test length(nodeiter(parsenewick("((,),(,));"))) == 7
    @test "MyLeaf" ∈ nodenameiter(parsenewick("((MyLeaf,),(,));"))
    tree = parsenewick("((MyLeaf:4.0,)Parent,(,));")
    branches = branchfilter(tree) do branch
        return src(branch) == "Parent" && dst(branch) == "MyLeaf"
    end
    @test length(branches) == 1
    @test getlength(first(branches)) ≈ 4.0
    @test_throws ErrorException parsenewick("((,),(,)));")
    @test_throws ErrorException parsenewick("((,),(,))")
    @test_throws ErrorException parsenewick("((,),(,);")
end

end
