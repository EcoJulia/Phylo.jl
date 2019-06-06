# Only run coverage from linux Julia 1.1 build on travis.
get(ENV, "TRAVIS_OS_NAME", "")       == "linux" || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "1.1"   || exit()

Pkg.add("Coverage")
using Coverage

cd(joinpath(dirname(@__FILE__), "..")) do
    processed = process_folder()
    Coveralls.submit(processed)
    Codecov.submit(processed)
end
