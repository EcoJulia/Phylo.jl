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

function nodedepths(tree::Phylo.AbstractTree)
    function finddepths!(clade::String, parentdepth::Float64 = 0.0)
        mydepth = parentdepth
        push!(names, clade)
        if hasinbound(tree, clade)
             mydepth += getlength(tree, getinbound(tree, clade))
        end
        depth[clade] = mydepth
        for ch in getchildren(tree, clade)
            finddepths!(ch, mydepth)
        end
    end

    depth = Dict{String, Float64}()
    names = String[]
    sizehint!(depth, nnodes(tree))
    sizehint!(names, nnodes(tree))
    root = getnodename(tree, getroot(tree))
    finddepths!(root)
    depth, names
end
