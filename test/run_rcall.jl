module ValidateRCall_ape
using Compat.Test
using Compat: @warn
using Phylo
using RCall

# Create a temporary directory to work in
libdir = mktempdir();

global skipR = false
if !rcopy(R"require(ape)")
    rcall(Symbol(".libPaths"), libdir);
    reval("install.packages(\"ape\", lib=\"$libdir\", " *
          "repos=\"http://cran.r-project.org\")");
    global skipR = !rcopy(R"require(ape, lib.loc=c(\"$libdir\", .libPaths()))") &&
        !mustCrossvalidate;
    skipR && @warn "ape R package not installed and would not install, " *
        "skipping R crossvalidation"
end

if !skipR
    # Run tree comparisons on increasing numbers of tips
    @testset "RCall - testing Phylo vs ape" begin
        @testset "Testing with R rtree($i)" for i in 5:5:50
            rt = rcall(:rtree, i)
            jt = rcopy(NamedTree, rt)
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

        @testset "Testing with julia rand(Nonultrametric($i))" for i in 5:5:50
            jt = rand(Nonultrametric(i))
            rt = RObject(jt)
            @test Set(getleafnames(jt)) ==
                Set(rcopy(rcall(Symbol("[["), rt, "tip.label")))
            jt2 = rcopy(NamedTree, rt)
            @test getleafnames(jt) == getleafnames(jt2)
            @test rcopy(rcall(Symbol("all.equal"), rt, RObject(jt2)))
        end

        @testset "Testing reading in newick trees from disk" begin
            jtree = open(parsenewick, Phylo.path("H1N1.newick"))
            rtree = rcall(Symbol("read.tree"), Phylo.path("H1N1.newick"))
            @test rcopy(rcall(Symbol("all.equal"), jtree, rtree))
        end

        @testset "Testing reading in nexus trees from disk" begin
            jts = open(parsenexus, Phylo.path("H1N1.trees"))
            rtree1 = R"read.nexus($(Phylo.path(\"H1N1.trees\")))$TREE1"
            jtree1 = jts["TREE1"]
            @test rcopy(rcall(Symbol("all.equal"), jtree1, rtree1))
            @test "H1N1_A_MIYAGI_3_2000" ∈ nodenameiter(jtree1)
            @test collect(keys(first(treeinfoiter(jts)))) == ["lnP"]
        end
    end
end
rm(libdir, force=true, recursive=true);

end
