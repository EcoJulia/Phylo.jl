using Compat
importall Phylo.API

"""
    BinaryNode{T}(AbstractVector{T}, AbstractVector{T}) <: AbstractNode

A node of strict binary phylogenetic tree
"""
type BinaryNode{T} <: AbstractNode
    inbound::Nullable{T}
    outbounds::Tuple{Nullable{T}, Nullable{T}}

    function (::Type{BinaryNode{T}}){T}(inbound::AbstractVector{T} = T[],
                                        outbounds::AbstractVector{T} = T[])
        length(inbound) <= 1 ||
            error("At most one inbound connection to BinaryNode")
        n_in = length(inbound) == 0 ? Nullable{T}() :
            Nullable(inbound[1])
        length(outbounds) <= 2 ||
            error("At most two outbound connections from BinaryNode")
        n_out = length(outbounds) == 0 ? (Nullable{T}(), Nullable{T}()) :
            (length(outbounds) == 1 ? (Nullable(outbounds[1]), Nullable{T}()) :
             (Nullable(outbounds[1]), Nullable(outbounds[2])))
        new{T}(n_in, n_out)
    end
end

function _hasinbound(node::BinaryNode)
    return !isnull(node.inbound)
end

function _outdegree(node::BinaryNode)
    return (isnull(node.outbounds[1]) ? 0 : 1) +
        (isnull(node.outbounds[2]) ? 0 : 1)
end

function _hasoutboundspace(node::BinaryNode)
    return _outdegree(node) < 2
end

function _getinbound(node::BinaryNode)
    _hasinbound(node) ||
        error("Node has no inbound connection")
    return get(node.inbound)
end

function _setinbound!{T}(node::BinaryNode{T}, inbound::T)
    !_hasinbound(node) ||
        error("BinaryNode already has an inbound connection")
    node.inbound = inbound
end

function _deleteinbound!{T}(node::BinaryNode{T}, inbound::T)
    _hasinbound(node) ||
        error("Node has no inbound connection")
    get(node.inbound) == inbound ||
        error("BinaryNode has no inbound connection from branch $inbound")
    node.inbound = Nullable{T}()
end

function _getoutbounds{T}(node::BinaryNode{T})
    return isnull(node.outbounds[1]) ?
        (isnull(node.outbounds[2]) ? T[] : [get(node.outbounds[2])]) :
        (isnull(node.outbounds[2]) ? [get(node.outbounds[1])] :
         [get(node.outbounds[1]), get(node.outbounds[2])])
end

function _addoutbound!{T}(node::BinaryNode{T}, outbound::T)
    isnull(node.outbounds[1]) ?
        node.outbounds = (Nullable(outbound), node.outbounds[2]) :
        (isnull(node.outbounds[2]) ?
         node.outbounds = (node.outbounds[1], Nullable(outbound)) :
         error("BinaryNode already has two outbound connections"))
end

function _deleteoutbound!{T}(node::BinaryNode{T}, outbound::T)
    !isnull(node.outbounds[1]) && get(node.outbounds[1]) == outbound ?
        node.outbounds = (node.outbounds[2], Nullable{T}()) :
        (!isnull(node.outbounds[2]) && get(node.outbounds[2]) == outbound ?
         node.outbounds = (node.outbounds[1], Nullable{T}()) :
         error("BinaryNode does not have outbound connection to branch $outbound"))
end
