module TestInference

# Environment variable to avoid boring R package builds
mustCrossvalidate = haskey(ENV, "JULIA_MUST_CROSSVALIDATE") &&
                    ENV["JULIA_MUST_CROSSVALIDATE"] == "1"

# Only run R on unix or when R is installed because JULIA_MUST_CROSSVALIDATE is set to 1
global skipR = !(mustCrossvalidate || Sys.isunix())
try
    skipR && error("Not on unix...")
    using RCall
    global skipR = false
catch
    global skipR = true
    if mustCrossvalidate
        error("R not installed, JULIA_MUST_CROSSVALIDATE set")
    else
        @warn "R not installed, skipping R cross-validation."
    end
end

!skipR && include("run_inference.jl")

end
