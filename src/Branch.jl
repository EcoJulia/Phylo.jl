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
                        len::Float64 = NaN) where {RT, NL}
        len >= 0.0 || isnan(len) ||
            error("Branch length must be positive or NaN (no recorded length), not $len")
        new{RT, NL}(name, source, destination, len)
    end
end

import Phylo.API: _src
_src(::AbstractTree, branch::Branch) = branch.source
import Phylo.API: _dst
_dst(::AbstractTree, branch::Branch) = branch.destination
import Phylo.API: _getlength
_getlength(::AbstractTree, branch::Branch) = branch.length
import Phylo.API: _getbranchname
_getbranchname(::AbstractTree, branch::Branch) = branch.name
import Phylo.API: _getbranchdata
_getbranchdata(::AbstractTree, ::Branch) = nothing

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
