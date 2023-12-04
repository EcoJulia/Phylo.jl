module ValidateRCall_wrapped
using Test
using Phylo
using RCall
using DataFrames

global skipR = !rcopy(R"require(ape)")

# Run tree comparisons on increasing numbers of tips
@testset "RCall - testing Phylo vs ape" begin
    @testset "For $TreeType" for TreeType in
        (skipR ? [] : [NamedTree, NamedBinaryTree,
                       BinaryTree{ManyRoots, DataFrame, Vector{Float64}},
                       Phylo.LTD{OneRoot, Float64}, Phylo.LTD{ManyRoots, Float64},
                       RootedTree, ManyRootTree])
                                              
        @testset "Testing with R rtree($i)" for i in 10:10:50
            rt = rcall(:rtree, i)
            jt = rcopy(TreeType, rt)
            jl = Set(getleafnames(jt))
            rl = Set(rcopy(rcall(Symbol("[["), rt, "tip.label")))
            @test jl == rl
            rt2 = RObject(jt)
            @test rcopy(rcall(Symbol("all.equal"), rt, rt2))
            @rput jt
            reval("jt2=jt");
            @rget jt2
            jt3=jt2;
            @rput jt3
            @rput rt
            @test rcopy(rcall(Symbol("all.equal"), rt, jt3))
        end

        @testset "Testing with julia rand(Nonultrametric($i))" for i in 10:10:50
            jt = rand(Nonultrametric{TreeType}(i))
            rt = RObject(jt)
            @test getleafnames(jt) == rcopy(rcall(Symbol("[["), rt, "tip.label"))
            jt2 = rcopy(TreeType, rt)
            @test Set(getleafnames(jt)) == Set(getleafnames(jt2))
            @test rcopy(rcall(Symbol("all.equal"), rt, RObject(jt2)))
        end

        @testset "Testing with julia rand(Ultrametric($i))" for i in 10:10:50
            jt = rand(Ultrametric{TreeType}(i))
            rt = RObject(jt)
            @test getleafnames(jt) == rcopy(rcall(Symbol("[["), rt, "tip.label"))
            jt2 = rcopy(TreeType, rt)
            @test Set(getleafnames(jt)) == Set(getleafnames(jt2))
            @test rcopy(rcall(Symbol("all.equal"), rt, RObject(jt2)))
        end

        @testset "Testing with julia rand(Ultrametric($i), [...])" for i in 10:10:50
            jt = rand(Ultrametric{TreeType}(i), ["one", "two"])
            rt = RObject(jt["one"])
            @test Set(getleafnames(jt)) ==
                Set(rcopy(rcall(Symbol("[["), rt, "tip.label")))
            jt2 = rcopy(TreeType, rt)
            @test Set(getleafnames(jt)) == Set(getleafnames(jt2))
            
            @test rcopy(rcall(Symbol("all.equal"), rt, RObject(jt2)))
        end

        @testset "Testing reading in newick trees from disk" begin
            if nodedatatype(TreeType) <: Dict
                jtree = open(io -> parsenewick(io, TreeType), Phylo.path("H1N1.newick"))
                rtree = rcall(Symbol("read.tree"), Phylo.path("H1N1.newick"))
                @test rcopy(rcall(Symbol("all.equal"), jtree, rtree))
            end
        end

        @testset "Testing reading in nexus trees from disk" begin
            if nodedatatype(TreeType) <: Dict
                jts = open(io -> parsenexus(io, TreeType), Phylo.path("H1N1.trees"))
                rtree1 = R"read.nexus($(Phylo.path(\"H1N1.trees\")))$TREE1"
                jtree1 = jts["TREE1"]
                @test rcopy(rcall(Symbol("all.equal"), jtree1, rtree1))
                @test "H1N1_A_MIYAGI_3_2000" ∈ nodenameiter(jtree1)
                @test collect(keys(gettreeinfo(jts)["TREE1"])) == ["lnP"]
            end
        end
    end
end

end
