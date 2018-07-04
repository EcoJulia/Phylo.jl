using Compat
using Compat: @info

# Identify files in test/ that are testing matching files in src/
#  - src/Source.jl will be matched by test/test_Source.jl

filebase = map(file -> replace(file, r"(.*).jl" => s"\1"),
                filter(file -> contains(file, r".*\.jl"), readdir("../src")))
testbase = map(file -> replace(file, r"test_(.*).jl" => s"\1"),
                filter(str -> contains(str, r"^test_.*\.jl$"), readdir()))

@info "Running tests for files:"
for t in testbase
    println("    = $t.jl")
end
println()

@info "Running tests..."
for t in testbase
    fn = "test_$t.jl"
    println("    * Testing $t.jl ...")
    include(fn)
    println()
end

# Identify tests with no matching file
superfluous = filter(f -> f ∉ filebase, testbase)
if length(superfluous) > 0
    @info "Potentially superfluous tests:"
    for f in superfluous
        println("    + $f.jl")
    end
    println()
end

# Identify files with no matching test
missing = filter(f -> f ∉ testbase, filebase)
if length(missing) > 0
    @info "Potentially missing tests:"
    for f in missing
        println("    - $f.jl")
    end
    println()
end
