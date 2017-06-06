using Compat
importall Phylo.API

"""
    Branch

    A directed branch connecting two AbstractNodes of phylogenetic tree
"""
type Branch{T} <: AbstractBranch
    source::T
    target::T
    length::Float64

    function (::Type{Branch{T}}){T}(source::T, target::T, length::Float64)
        length >= 0.0 || isnan(length) ||
            error("Branch length must be positive or NaN (no recorded length)")
        new{T}(source, target, length)
    end
end

Branch{T}(source::T, target::T, length::Float64) =
    Branch{T}(source, target, length)

const SimpleBranch = Branch{Int}

_getsource(branch::Branch) = branch.source
_gettarget(branch::Branch) = branch.target
_setsource!{T}(branch::Branch{T}, source::T) = branch.source = source
_settarget!{T}(branch::Branch{T}, target::T) = branch.target = target
_getlength(branch::Branch) = branch.length

function checkbranch(id::Int, branch::Branch, tree::AbstractTree)
    return id > 0 &&
        getsource(branch) != gettarget(branch) &&
        !haskey(getbranches(tree), id) &&
        haskey(getnodes(tree), getsource(branch)) &&
        haskey(getnodes(tree), gettarget(branch)) &&
        !hasinbound(getnodes(tree)[gettarget(branch)]) &&
        outboundspace(getnodes(tree)[getsource(branch)])
end
