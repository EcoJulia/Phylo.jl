# SPDX-License-Identifier: BSD-2-Clause

module TestPhylo
using Test
using Phylo
using Git
using Logging
using Pkg

function is_repo_clean(repo_path)
    # Get the status of the repository
    statuses = readlines(`$(Git.git()) status -s $repo_path`)
    statuses = filter((!) ∘ contains("??"), statuses)

    is_clean = isempty(statuses)
    is_clean || @error "\n" * join(statuses, "\n")

    return is_clean
end

if !haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows"
    Pkg.develop(url = "https://github.com/richardreeve/ResearchSoftwareMetadata.jl.git")
    using ResearchSoftwareMetadata

    @testset "RSMD" begin
        git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
        @test isnothing(ResearchSoftwareMetadata.crosswalk())
        global_logger(SimpleLogger(stderr, Logging.Warn))
        @test_nowarn ResearchSoftwareMetadata.crosswalk()
        @test is_repo_clean(git_dir)
    end
end

@testset "Deprecations" begin
    t = NamedTree()
    @test_deprecated name = addnode!(t)
    @test_deprecated addnode!(t, "Second")
    @test_deprecated addbranch!(t, name, "Second")
    @test_deprecated branch!(t, name)
    @test nleaves(t) == 2
    @test nroots(t) == 1
    @test nbranches(t) == 2
    @test_deprecated treeiter(t)
    @test_deprecated treenameiter(t)
end

end
