import Phylo.API: _src, _dst, _setsrc!, _setdst!, _getlength

"""
    Branch

    A directed branch connecting two AbstractNodes of phylogenetic tree
"""
mutable struct Branch{T} <: AbstractBranch
    source::T
    destination::T
    length::Float64

    function Branch(source::T, destination::T, length::Float64) where T
        length >= 0.0 || isnan(length) ||
            error("Branch length must be positive or NaN (no recorded length)")
        new{T}(source, destination, length)
    end
end

const SimpleBranch = Branch{Int}

_src(branch::Branch) = branch.source
_dst(branch::Branch) = branch.destination
_setsrc!(branch::Branch{T}, source::T) where T = branch.source = source
_setdst!(branch::Branch{T}, destination::T) where T = branch.destination = destination
_getlength(branch::Branch) = branch.length

function checkbranch(id::Int, branch::Branch, tree::AbstractTree)
    return id > 0 &&
        src(branch) != dst(branch) &&
        !haskey(getbranches(tree), id) &&
        haskey(getnodes(tree), src(branch)) &&
        haskey(getnodes(tree), dst(branch)) &&
        !hasinbound(getnodes(tree)[dst(branch)]) &&
        outboundspace(getnodes(tree)[src(branch)])
end
