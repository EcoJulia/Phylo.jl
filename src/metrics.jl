function mrca(tree, target)
    ancestors = getancestors(tree, first(target))
    ranks = Dict(j => i for (i,j) in enumerate(ancestors))
    checked = Set(ancestors)
    oldest = 1
    for species in target
        while !(species ∈ checked)
            push!(checked, species)
            species = getparent(tree, species)
        end
        oldest = max(oldest, get(ranks, species, 0))
    end
    ancestors[oldest]
end
