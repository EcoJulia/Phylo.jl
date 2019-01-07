import Phylo.API: _src, _dst, _getlength, _getbranchname

"""
    Branch

    A branch connecting two AbstractNodes of a phylogenetic tree
"""
mutable struct Branch{RT, NL} <: AbstractBranch{RT, NL}
    name::Int
    source::NL
    destination::NL
    length::Float64

    function Branch{RT}(name::Int, source::NL, destination::NL,
                        length::Float64 = NaN) where {RT, NL}
        length >= 0.0 || isnan(length) ||
            error("Branch length must be positive or NaN (no recorded length), not $length")
        new{RT, NL}(name, source, destination, length)
    end
end

_src(::AbstractTree, branch::Branch) = branch.source
_dst(::AbstractTree, branch::Branch) = branch.destination
_getlength(::AbstractTree, branch::Branch) = branch.length
_getbranchname(::AbstractTree, branch::Branch) = branch.name

function checkbranch(id::Int, branch::Branch{<: Rooted}, tree::AbstractTree)
    return id > 0 &&
        hasbranch(tree, id) &&
        src(tree, branch) != dst(tree, branch) &&
        hasnode(tree, src(tree, branch)) &&
        hasnode(tree, dst(tree, branch)) &&
        indegree(tree, getnode(tree, dst(tree, branch))) > 0 &&
        outdegree(tree, src(tree, branch)) > 0 &&
        getbranch(tree, getinbound(tree, dst(tree, branch))) === branch
end
