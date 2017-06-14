module TestRCall_ape
using Base.Test

using Phylo
Rinstalled = false

module addmacros
macro rput(x) end
macro rget(x) end

export @rput, @rget
end

# Only run R on unix (macs and linux)
skipR = !is_unix()

# Environment variable to avoid boring R package builds
skipRinstall = haskey(ENV, "SKIP_R_INSTALL") && ENV["SKIP_R_INSTALL"] == "1"

try
    skipR && error("Skipping R testing...")
    using RCall
    include(joinpath(dirname(dirname(@__FILE__)), "src", "rcall.jl"))
    Rinstalled = true
catch
    warn("R not installed, skipping RCall testing")
    using TestRCall_ape.addmacros
end

if Rinstalled
    # Create a temporary directory to work in
    libdir = mktempdir();

    # Only run R on macs
    if !skipR
        # Skip the (slow!) R package installation step
        if skipRinstall
            reval("library(ape)");
        else
            rcall(Symbol(".libPaths"), libdir);
            reval("install.packages(c(\"devtools\", \"ggplot2\", \"ape\", \"plyr\", \"phangorn\", \"tidyr\", \"tibble\", \"phytools\", \"reshape2\", \"ggthemes\"), lib=\"$libdir\", repos=\"http://cran.r-project.org\", type=\"binary\")");
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
        end
    end
    rm(libdir, force=true, recursive=true);
end

end
