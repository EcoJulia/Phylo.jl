using Compat

"""
    BinaryNode{B}(AbstractVector{B}, AbstractVector{B}) <: AbstractNode

A node of strict binary phylogenetic tree
"""
mutable struct BinaryNode{RT <: Rooted, NL,
                          B <: AbstractBranch{RT, NL}} <: AbstractNode{RT, NL}
    name::NL
    inbound::Union{B, Nothing}
    outbounds::Tuple{Union{B, Nothing}, Union{B, Nothing}}

    function BinaryNode{RT, NL, B}(name::NL,
                                   inbound::AbstractVector{B} = B[],
                                   outbounds::AbstractVector{B} = B[]) where
        {RT, NL, B <: AbstractBranch{RT, NL}}
        length(inbound) <= 1 ||
            error("At most one inbound connection to BinaryNode")
        n_in = length(inbound) == 0 ? nothing : inbound[1]
        length(outbounds) <= 2 ||
            error("At most two outbound connections from BinaryNode")
        n_out = length(outbounds) == 0 ? (nothing, nothing) :
            (length(outbounds) == 1 ? (outbounds[1], nothing) :
             (outbounds[1], outbounds[2]))
        new{RT, NL, B}(name, n_in, n_out)
    end
end

import Phylo.API._hasinbound
function _hasinbound(::AbstractTree, node::BinaryNode)
    return node.inbound != nothing
end

import Phylo.API._outdegree
function _outdegree(::AbstractTree{OneTree, RT, NL},
                    node::BinaryNode{RT, NL}) where {RT <: Rooted, NL}
    return (node.outbounds[1] === nothing ? 0 : 1) +
        (node.outbounds[2] === nothing ? 0 : 1)
end

import Phylo.API._hasoutboundspace
_hasoutboundspace(tree::AbstractTree{OneTree, RT, NL, N, B},
                  node::N) where {RT <: Rooted, NL,
                                  N <: BinaryNode{RT, NL}, B} =
    _outdegree(tree, node) < 2
_hasoutboundspace(tree::AbstractTree{OneTree, RT, NL, N, B},
                  node::NL) where {RT <: Rooted, NL,
                                   N <: BinaryNode{RT, NL}, B} =
    _outdegree(tree, _getnode(tree, node)) < 2


import Phylo.API._getinbound
function _getinbound(tree::AbstractTree, node::BinaryNode)
    _hasinbound(tree, node) ||
        error("Node has no inbound connection")
    return node.inbound
end

import Phylo.API._addinbound!
function _addinbound!(tree::AbstractTree, node::BinaryNode{RT, NL, B},
                      inbound::B) where {RT, NL, B <: AbstractBranch}
    _hasinbound(tree, node) &&
        error("BinaryNode already has an inbound connection")
    node.inbound = inbound
end

import Phylo.API._removeinbound!
function _removeinbound!(tree::AbstractTree,
                         node::BinaryNode{RT, NL, B},
                         inbound::B) where
    {RT, NL, B <: AbstractBranch}
    _hasinbound(tree, node) || error("Node has no inbound connection")
    node.inbound == inbound ||
        error("BinaryNode has no inbound connection from branch $inbound")
    node.inbound = nothing
end

import Phylo.API._getoutbounds
function _getoutbounds(::AbstractTree, node::BinaryNode)
    return (node.outbounds[1] === nothing ?
        (node.outbounds[2] === nothing ? T[] : [node.outbounds[2]]) :
        (node.outbounds[2] === nothing ? [node.outbounds[1]] :
         [node.outbounds[1], node.outbounds[2]]))
end

import Phylo.API._addoutbound!
function _addoutbound!(::AbstractTree, node::BinaryNode{RT, NL, B},
                       outbound::B) where {RT, NL, B <: AbstractBranch}
    node.outbounds[1] === nothing ?
        node.outbounds = (outbound, node.outbounds[2]) :
        (node.outbounds[2] === nothing ?
         node.outbounds = (node.outbounds[1], outbound) :
         error("BinaryNode already has two outbound connections"))
end

import Phylo.API._removeoutbound!
function _removeoutbound!(::AbstractTree, node::BinaryNode{RT, NL, B},
                          outbound::B) where
    {RT, NL, B <: AbstractBranch}
    node.outbounds[1] == outbound ?
        node.outbounds = (node.outbounds[2], nothing) :
        (node.outbounds[2] == outbound ?
         node.outbounds = (node.outbounds[1], nothing) :
         error("BinaryNode does not have outbound connection to branch " *
               "$outbound"))
end

import Phylo.API._getnodename
_getnodename(::AbstractTree, node::BinaryNode) = node.name

"""
    Node{RT, NL, T}(AbstractVector{T}, AbstractVector{T}) <: AbstractNode

A node of potentially polytomous phylogenetic tree
"""
mutable struct Node{RT <: Rooted, NL, B <: AbstractBranch{RT, NL}} <:
    AbstractNode{RT, NL}
    name::NL
    inbound::Union{B, Nothing}
    outbounds::Vector{B}

    function Node{RT, NL, B}(name::NL,
                             inbound::AbstractVector{B} = B[],
                             outbounds::AbstractVector{B} = B[]) where
        {RT, NL, B <: AbstractBranch{RT, NL}}
        length(inbound) <= 1 ||
            error("At most one inbound connection to Node")
        n_in = length(inbound) == 0 ? nothing : inbound[1]
        n_out = outbounds
        new{RT, NL, B}(name, n_in, n_out)
    end
end

function _hasinbound(::AbstractTree, node::Node)
    return node.inbound != nothing
end

function _outdegree(::AbstractTree{OneTree, RT, NL},
                    node::Node{RT, NL}) where {RT <: Rooted, NL}
    return length(node.outbounds)
end

_hasoutboundspace(tree::AbstractTree{OneTree, RT, NL, N, B},
                  node::N) where {RT <: Rooted, NL,
                                  N <: Node{RT, NL}, B} = true
_hasoutboundspace(tree::AbstractTree{OneTree, RT, NL, N, B},
                  node::NL) where {RT <: Rooted, NL,
                                   N <: Node{RT, NL}, B} = true

import Phylo.API._getinbound
function _getinbound(tree::AbstractTree, node::Node)
    _hasinbound(tree, node) ||
        error("Node has no inbound connection")
    return node.inbound
end

import Phylo.API._addinbound!
function _addinbound!(tree::AbstractTree, node::Node{RT, NL, B},
                      inbound::B) where {RT, NL, B <: AbstractBranch}
    _hasinbound(tree, node) &&
        error("Node already has an inbound connection")
    node.inbound = inbound
end

import Phylo.API._removeinbound!
function _removeinbound!(tree::AbstractTree,
                         node::Node{RT, NL, B},
                         inbound::B) where
    {RT, NL, B <: AbstractBranch}
    _hasinbound(tree, node) || error("Node has no inbound connection")
    node.inbound == inbound ||
        error("Node has no inbound connection from branch $inbound")
    node.inbound = nothing
end

import Phylo.API._getoutbounds
function _getoutbounds(::AbstractTree, node::Node)
    return node.outbounds
end

import Phylo.API._addoutbound!
function _addoutbound!(::AbstractTree, node::Node{RT, NL, B},
                       outbound::B) where {RT, NL, B <: AbstractBranch}
    push!(node.outbounds, outbound)
end

import Phylo.API._removeoutbound!
function _removeoutbound!(::AbstractTree, node::Node{RT, NL, B},
                          outbound::B) where {RT, NL,
                                              B <: AbstractBranch}
    outbound ∈ node.outbounds ? filter!(n -> n != outbound, node.outbounds) :
         error("Node does not have outbound connection to branch $outbound")
end

import Phylo.API._getnodename
_getnodename(::AbstractTree, node::Node) = node.name
