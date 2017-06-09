module TestRCall_ape
using Base.Test

using Phylo
Rinstalled = false
try
    using RCall
    Rinstalled = true
    include("../src/rcall.jl")
catch
    warn("R not installed, skipping RCall testing")
end

if Rinstalled
    # Create a temporary directory to work in
    libdir = mktempdir();

    # Only run R on macs
    if is_apple()
        # Environment variable to avoid boring R package builds
        skipRinstall = haskey(ENV, "SKIP_R_INSTALL") && ENV["SKIP_R_INSTALL"] == "1"
        # Skip the (slow!) R package installation step
        if skipRinstall
            reval("library(ape)");
        else
            rcall(Symbol(".libPaths"), libdir);
            rcall(Symbol("install.packages"),
                  ["devtools", "ggplot2", "ape",
                   "phangorn", "tidyr", "tibble", "phytools",
                   "reshape2", "ggthemes"],
                  lib=libdir, repos="http://cran.r-project.org");
            reval("library(ape, lib.loc=c(\"$libdir\", .libPaths()))");
        end

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
        end
    end
    rm(libdir, force=true, recursive=true);
end

end
