# SPDX-License-Identifier: BSD-2-Clause

module TestPhylo
using Test
using Phylo
using Git
using Logging
using Pkg

function is_repo_clean(repo_path; ignore_untracked = true)
    # Get the status of the repository
    statuses = readlines(`$(Git.git()) status -s $repo_path`)

    if ignore_untracked
        # Repo must be clean except for untracked files
        statuses = filter((!) ∘ contains("??"), statuses)
    end

    is_clean = isempty(statuses)

    # If not clean then report on dirty files
    is_clean || @error "\n" * join(statuses, "\n")

    return is_clean
end

# Metadata crosswalk testing only works on Julia v1.8 and after due to Project.toml changes
# Also does not currently work on Windows runners on GitHub due to file writing issues
if VERSION ≥ VersionNumber("1.8.0") && (!haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows")
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
