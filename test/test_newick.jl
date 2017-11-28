module TestNewick

using Phylo
using Compat.Test

@testset "A few simple trees" begin
    @test length(nodeiter(parsenewick("((,),(,));"))) == 7
    @test ["where", "when it's good", "Not mine", "MyLeaf"] ⊆
        nodenameiter(parsenewick("""((MyLeaf,"when it's good"),
                                     ('Not mine',where));"""))
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
