using Compat
using Compat.Test
using Compat: @info

# Identify files in test/ that are testing matching files in src/
#  - src/Source.jl will be matched by test/test_Source.jl

filebase = map(file -> replace(file, r"(.*).jl$" => s"\1"),
                filter(file -> occursin(r".*\.jl$", file),
                       readdir("../src")))
testbase = map(file -> replace(file, r"test_(.*).jl$" => s"\1"),
                filter(str -> occursin(r"^test_.*\.jl$", str), readdir()))

println()
@info "Running tests for files:"
for t in testbase
    println("    = $t.jl")
end
println()

@testset "* Testing $t.jl" for t in testbase
    include("test_$t.jl")
end

# Identify tests with no matching file
superfluous = filter(f -> f ∉ filebase, testbase)
if length(superfluous) > 0
    println()
    @info "Potentially superfluous tests:"
    for f in superfluous
        println("    + $f.jl")
    end
    println()
end

# Identify files with no matching test
missing = filter(f -> f ∉ testbase, filebase)
if length(missing) > 0
    println()
    @info "Potentially missing tests:"
    for f in missing
        println("    - $f.jl")
    end
    println()
end

# Identify files that are cross-validating results against other packages
# test/pkg_Package.jl should validate results against the Package package

pkgbase = map(file -> replace(file, r"pkg_(.*).jl$" => s"\1"),
                   filter(str -> occursin(r"^pkg_.*\.jl$", str),
                          readdir()))

if length(pkgbase) > 0
    @info "Running cross-validation against:"
    for p in pkgbase
        println("    = $p")
    end
    println()

    @testset " * Validating against $p" for p in pkgbase
        include("pkg_$p.jl")
    end
end
