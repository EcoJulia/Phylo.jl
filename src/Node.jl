using Compat

"""
    BinaryNode{T}(AbstractVector{T}, AbstractVector{T}) <: AbstractNode

A node of strict binary phylogenetic tree
"""
mutable struct BinaryNode{T} <: AbstractNode
    inbound::Union{T, Nothing}
    outbounds::Tuple{Union{T, Nothing}, Union{T, Nothing}}

    function BinaryNode{T}(inbound::AbstractVector{T} = T[],
                           outbounds::AbstractVector{T} = T[]) where T
        length(inbound) <= 1 ||
            error("At most one inbound connection to BinaryNode")
        n_in = length(inbound) == 0 ? nothing : inbound[1]
        length(outbounds) <= 2 ||
            error("At most two outbound connections from BinaryNode")
        n_out = length(outbounds) == 0 ? (nothing, nothing) :
            (length(outbounds) == 1 ? (outbounds[1], nothing) :
             (outbounds[1], outbounds[2]))
        new{T}(n_in, n_out)
    end
end

import Phylo.API._hasinbound
function _hasinbound(node::BinaryNode)
    return node.inbound != nothing
end


import Phylo.API._outdegree
function _outdegree(node::BinaryNode)
    return (node.outbounds[1] === nothing ? 0 : 1) +
        (node.outbounds[2] === nothing ? 0 : 1)
end


import Phylo.API._hasoutboundspace
function _hasoutboundspace(node::BinaryNode)
    return _outdegree(node) < 2
end

import Phylo.API._getinbound
function _getinbound(node::BinaryNode)
    _hasinbound(node) ||
        error("Node has no inbound connection")
    return node.inbound
end

import Phylo.API._setinbound!
function _setinbound!(node::BinaryNode{T}, inbound::T) where T
    !_hasinbound(node) ||
        error("BinaryNode already has an inbound connection")
    node.inbound = inbound
end

import Phylo.API._deleteinbound!
function _deleteinbound!(node::BinaryNode{T}, inbound::T) where T
    _hasinbound(node) ||
        error("Node has no inbound connection")
    node.inbound == inbound ||
        error("BinaryNode has no inbound connection from branch $inbound")
    node.inbound = nothing
end

import Phylo.API._getoutbounds
function _getoutbounds(node::BinaryNode{T}) where T
    return node.outbounds[1] === nothing ?
        (node.outbounds[2] === nothing ? T[] : [node.outbounds[2]]) :
        (node.outbounds[2] === nothing ? [node.outbounds[1]] :
         [node.outbounds[1], node.outbounds[2]])
end

import Phylo.API._addoutbound!
function _addoutbound!(node::BinaryNode{T}, outbound::T) where T
    node.outbounds[1] === nothing ?
        node.outbounds = (outbound, node.outbounds[2]) :
        (node.outbounds[2] === nothing ?
         node.outbounds = (node.outbounds[1], outbound) :
         error("BinaryNode already has two outbound connections"))
end

import Phylo.API._deleteoutbound!
function _deleteoutbound!(node::BinaryNode{T}, outbound::T) where T
    node.outbounds[1] == outbound ?
        node.outbounds = (node.outbounds[2], nothing) :
        (node.outbounds[2] == outbound ?
         node.outbounds = (node.outbounds[1], nothing) :
         error("BinaryNode does not have outbound connection to branch $outbound"))
end


"""
    Node{T}(AbstractVector{T}, AbstractVector{T}) <: AbstractNode

A node of potentially polytomous phylogenetic tree
"""
mutable struct Node{T} <: AbstractNode
    inbound::Union{T, Nothing}
    outbounds::Vector{T}

    function Node{T}(inbound::AbstractVector{T} = T[],
                           outbounds::AbstractVector{T} = T[]) where T
        length(inbound) <= 1 ||
            error("At most one inbound connection to Node")
        n_in = length(inbound) == 0 ? nothing : inbound[1]
        n_out = outbounds
        new{T}(n_in, n_out)
    end
end

import Phylo.API._hasinbound
function _hasinbound(node::Node)
    return node.inbound != nothing
end


import Phylo.API._outdegree
function _outdegree(node::Node)
    return length(node.outbounds)
end


import Phylo.API._hasoutboundspace
function _hasoutboundspace(node::Node)
    return true
end

import Phylo.API._getinbound
function _getinbound(node::Node)
    _hasinbound(node) ||
        error("Node has no inbound connection")
    return node.inbound
end

import Phylo.API._setinbound!
function _setinbound!(node::Node{T}, inbound::T) where T
    !_hasinbound(node) ||
        error("Node already has an inbound connection")
    node.inbound = inbound
end

import Phylo.API._deleteinbound!
function _deleteinbound!(node::Node{T}, inbound::T) where T
    _hasinbound(node) ||
        error("Node has no inbound connection")
    node.inbound == inbound ||
        error("Node has no inbound connection from branch $inbound")
    node.inbound = nothing
end

import Phylo.API._getoutbounds
function _getoutbounds(node::Node{T}) where T
    return node.outbounds
end

import Phylo.API._addoutbound!
function _addoutbound!(node::Node{T}, outbound::T) where T
    push!(node.outbounds, outbound)
end

import Phylo.API._deleteoutbound!
function _deleteoutbound!(node::Node{T}, outbound::T) where T
    outbound ∈ node.outbounds ? filter!(n -> n != outbound, node.outbounds) :
         error("Node does not have outbound connection to branch $outbound")
end
