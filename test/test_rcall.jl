module ValidateRCall_rcall
macro R_str(script)
    script
end
macro rput(script)
    script
end
macro rget(script)
    script
end
export @R_str, @rput, @rget
end

module ValidateRCall_ape
using Compat.Test
using Compat: @warn

using Phylo

# Environment variable to avoid boring R package builds
mustCrossvalidate = haskey(ENV, "JULIA_MUST_CROSSVALIDATE") && ENV["JULIA_MUST_CROSSVALIDATE"] == "1"

# Only run R on unix or when R is installed because JULIA_MUST_CROSSVALIDATE is set to 1
skipR = !mustCrossvalidate && !is_unix()
Rinstalled = false
try
    skipR && error("Skipping R testing...")
    using RCall
    include(joinpath(dirname(dirname(@__FILE__)), "src", "rcall.jl"))
    Rinstalled = true
catch
    if mustCrossvalidate
        error("R not installed, but JULIA_MUST_CROSSVALIDATE is set")
    else
        @warn "R or appropriate Phylo package not installed, skipping R cross-validation."
    end
    using ValidateRCall_rcall
end

if Rinstalled
    # Create a temporary directory to work in
    libdir = mktempdir();

    if !skipR && !rcopy(R"require(ape)")
        rcall(Symbol(".libPaths"), libdir);
        reval("install.packages(\"ape\", lib=\"$libdir\", " *
              "repos=\"http://cran.r-project.org\")");
        skipR = !rcopy(R"require(ape, lib.loc=c(\"$libdir\", .libPaths()))") &&
            !mustCrossvalidate;
        skipR && @warn "ape R package not installed and would not install, " *
            "skipping R crossvalidation"
    end

    if !skipR
        # Run tree comparisons on increasing numbers of tips
        @testset "RCall - testing Phylo vs ape" begin
            @testset "Testing with R rtree($i)" for i in 5:5:50
                rt = rcall(:rtree, i)
                jt = NamedTree(rt)
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
                jt2 = NamedTree(rt)
                @test getleafnames(jt) == getleafnames(jt2)
                @test rcopy(rcall(Symbol("all.equal"), rt, RObject(jt2)))
            end

            @testset "Testing reading in newick trees from disk" begin
                jtree = open(parsenewick, "h1n1.trees")
                rtree = rcall(Symbol("read.tree"), "h1n1.trees")
                @test rcopy(rcall(Symbol("all.equal"), jtree, rtree))
            end
        end
    end
    rm(libdir, force=true, recursive=true);
end

end
