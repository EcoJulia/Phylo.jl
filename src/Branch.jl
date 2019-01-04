import Phylo.API: _src, _dst, _getlength

"""
    Branch

    A branch connecting two AbstractNodes of a phylogenetic tree
"""
mutable struct Branch{RT, NL} <: AbstractBranch{RT, NL}
    source::NL
    destination::NL
    length::Float64

    function Branch{RT}(source::NL, destination::NL,
                        length::Float64 = NaN) where {RT, NL}
        length >= 0.0 || isnan(length) ||
            error("Branch length must be positive or NaN (no recorded length)")
        new{RT, NL}(source, destination, length)
    end
end

_src(::AbstractTree, branch::Branch) = branch.source
_dst(::AbstractTree, branch::Branch) = branch.destination
_getlength(::AbstractTree, branch::Branch) = branch.length

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
