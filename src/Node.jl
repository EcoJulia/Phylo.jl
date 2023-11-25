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

import Phylo.API: _prefernodeobjects
_prefernodeobjects(::Type{<:BinaryNode}) = false

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

_prefernodeobjects(::Type{<:Node}) = false
