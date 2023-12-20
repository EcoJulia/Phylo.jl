using Unitful

"""
    struct RecursiveElt <: AbstractElt

A type for branches or nodes in a RecursiveTree, allowing navigation of the tree without
using the tree object itself.
"""
Base.@kwdef mutable struct RecursiveElt{RT, NL, MyType <: AbstractElt{RT, NL}, MyData,
                                        TheirType <: AbstractElt{RT, NL}, TheirData,
                                        BT <: BranchingType, LenUnits <: Number} <:
    AbstractElt{RT, NL}

    name::Union{NL, Nothing} = nothing
    id::Union{Int, Missing} = missing
    in::Union{RecursiveElt{RT, NL, TheirType, TheirData, MyType, MyData,
                           BT, LenUnits}, Nothing} = nothing
    conns::Vector{RecursiveElt{RT, NL, TheirType, TheirData, MyType, MyData,
                               BT, LenUnits}} =
        RecursiveElt{RT, NL, TheirType, TheirData, MyType, MyData, BT, LenUnits}[]
    data::MyData = _emptydata(MyData)
    length::Union{LenUnits, Missing} = missing
end

const RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits} =
    RecursiveElt{RT, NL, AbstractNode{RT, NL}, NodeData,
                 AbstractBranch{RT, NL}, BranchData, BT, LenUnits}

const RecursiveBranch{RT, NL, NodeData, BranchData, BT, LenUnits} =
    RecursiveElt{RT, NL, AbstractBranch{RT, NL}, BranchData,
                 AbstractNode{RT, NL}, NodeData, BT, LenUnits}    

RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits}(name::NL, data::NodeData = nothing) where
    {RT, NL, NodeData, BranchData, BT, LenUnits} =
    RecursiveNode{RT, NL, NodeData, BranchData,
                  BT, LenUnits}(name, nothing,
                                RecursiveBranch{RT, NL, NodeData,
                                                BranchData, LenUnits}[],
                                data, missing)

"""
    struct RecoursiveTree <: AbstractTree

A phylogenetic tree type containing RecursiveElts as both nodes and branches,
allowing navigation of the tree  using only the node and branch elements.
"""
mutable struct RecursiveTree{RT, NL, NodeData, BranchData, BT <: BranchingType, LenUnits, TD} <:
    AbstractTree{OneTree, RT, NL,
                 RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits},
                 RecursiveBranch{RT, NL, NodeData, BranchData, BT, LenUnits}}

    name::String
    nodedict::Dict{NL, Int}
    roots::Vector{RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits}}
    nodes::Vector{Union{RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits},
                        Missing}}
    branches::Vector{Union{RecursiveBranch{RT, NL, NodeData, BranchData, BT, LenUnits},
                           Missing}}
    data::Dict{String, Any}
    tipdata::TD
    rootheight::Union{LenUnits, Missing}
    isvalid::Union{Bool, Missing}
    cache::Dict{TraversalOrder,
                Vector{RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits}}}

    function RecursiveTree{RT, NL, NodeData, BranchData,
                           BT, LenUnits, TD}(tipnames::AbstractVector{NL} = NL[];
                                             name::String = TREENAME,
                                             tipdata::TD = _emptydata(TD),
                                             rootheight::Union{LenUnits, Missing} = missing,
                                             validate::Bool = false) where
        {RT, NL, NodeData, BranchData, BT, LenUnits, TD}
        NT = RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits}
        BrT = RecursiveBranch{RT, NL, NodeData, BranchData, BT, LenUnits}
        tree = new{RT, NL, NodeData, BranchData,
                   BT, LenUnits, TD}(name, Dict{NL, NT}(), NT[],
                                     Union{NT, Missing}[], Union{BrT, Missing}[],
                                     Dict{String, Any}(),
                                     tipdata, rootheight, missing,
                                     Dict{TraversalOrder, Vector{NT}}())
        
        if !isempty(tipnames)
            createnodes!(tree, tipnames)
        elseif !isnothing(tree.tipdata) && !isempty(tree.tipdata)
            createnodes!(tree, unique(keys(tree.tipdata)))
        end

        if validate
            validate!(tree)
        else
            tree.isvalid = missing
        end

        return tree
    end
end

function RecursiveTree{RT, NL, NodeData, BranchData, BT, LenUnits, TD}(leafinfos::TD) where
    {RT, NL, NodeData, BranchData, BT, LenUnits, TD}

    leafnames = unique(info[1] for info in getiterator(leafinfos))
    return RecursiveTree{RT, NL, NodeData, BranchData, BT, LenUnits, TD}(leafnames; tipdata = leafinfos)
end

import Phylo.API: _prefernodeobjects, _preferbranchobjects
_prefernodeobjects(::Type{<:RecursiveTree}) = true
_prefernodeobjects(::Type{<:RecursiveNode}) = true
_prefernodeobjects(::Type{<:RecursiveElt}) = true
_preferbranchobjects(::Type{<:RecursiveTree}) = true
_preferbranchobjects(::Type{<:RecursiveBranch}) = true
_preferbranchobjects(::Type{<:RecursiveElt}) = true

import Phylo.API: _validate!
function _validate!(tree::RecursiveTree{RT, NL, NodeData, BranchData,
                                        BT, LenUnits, TD}) where
    {RT, NL, NodeData, BranchData, BT, LenUnits, TD}

    tree.isvalid = true

    if !_matchlabels(NL, _tiplabeltype(TD))
        tree.isvalid = false
        error("Tree $(_gettreename(tree)) has inconsistent node and tip label types $NL and $(_tiplabeltype(TD))")
    end

    if !isnothing(tree.tipdata) && !isempty(tree.tipdata)
        tree.isvalid &= (Set(keys(tree.tipdata)) == Set(getleafnames(tree)))
    end

    nr = nroots(tree)
    if RT ≡ OneRoot
        if nr ≠ 1
            @warn "Wrong number of roots for $RT tree ($nr)"
            tree.isvalid = false
        end
    elseif RT ≡ ManyRoots
        if nr < 1
            @warn "Wrong number of roots for $RT tree ($nr)"
            tree.isvalid = false
        end
    end

    if length(tree.nodedict) ≠ _nnodes(tree)
        @warn "Number of nodes in node lookup is inconsistent with " *
            "node vector for tree ($(length(tree.nodedict)) ≠ $(_nnodes(tree))"
        tree.isvalid = false
    else
        for (i, node) in enumerate(tree.nodes)
            if !ismissing(node)
                if i ≠ node.id || i ≠ tree.nodedict[node.name]
                    @warn "Mismatch between node $(node.name) id $(node.id) and " *
                        "references in tree: $i and $(get(tree.nodedict, node.name, missing))"
                    tree.isvalid = false
                else
                    for branch in node.conns
                        if node ≢ branch.in && node ∉ branch.conns
                            tree.isvalid = false
                            @warn "Mismatch between node name $(node.name) and its connections"
                        end
                    end
                    if RT <: Rooted
                        if _hasinbound(tree, node)
                            branch = node.in
                            if node ≢ branch.in && node ∉ branch.conns
                                tree.isvalid = false
                                @warn "Mismatch between node name $(node.name) and its connections"
                            end
                        end
                        if BT ≡ BinaryBranching && length(getoutbounds(tree, node)) > 2
                            tree.isvalid = false
                            @warn "Node name $(node.name) has too many outbound connections"
                        end
                    else # Unrooted
                        if BT ≡ BinaryBranching && length(getoutbounds(tree, node)) > 3
                            tree.isvalid = false
                            @warn "Node name $(node.name) has too many outbound connections"
                        end
                    end
                end
            end
        end
    end

    for key in keys(tree.cache)
        if length(tree.cache[key]) ≠ _nnodes(tree)
            @warn "Number of nodes in tree cache with key $key is inconsistent " *
                "with tree ($(length(tree.cache[key])) ≠ $(_nnodes(tree)))"
            tree.isvalid = false
        end
    end

    for (i, branch) in enumerate(tree.branches)
        if !ismissing(branch)
            if i ≠ branch.id
                @warn "Mismatch between branch id $(branch.id) and reference in tree $i"
                tree.isvalid = false
            else
                for node in branch.conns
                    if branch ≢ node.in && branch ∉ node.conns
                        tree.isvalid = false
                        @warn "Mismatch between branch id $(branch.id) and its connections"
                    end
                end
                if RT <: Rooted
                    node = branch.in
                    if branch ≢ node.in && branch ∉ node.conns
                        tree.isvalid = false
                        @warn "Mismatch between branch id $(branch.id) and its connections"
                    end
                end
            end
        end
    end

    return tree.isvalid
end

import Phylo.API: _invalidate!
function _invalidate!(tree::RecursiveTree, state = missing)
    empty!(tree.cache)
    tree.isvalid = state
end

# Type aliases

const ReB{RT, BT, LenUnits} = RecursiveBranch{RT, String, Dict{String, Any}, Dict{String, Any}, BT, LenUnits}
const ReN{RT, BT, LenUnits} = RecursiveNode{RT, String, Dict{String, Any}, Dict{String, Any}, BT, LenUnits}
const ReT{RT, TD, BT, LenUnits} = RecursiveTree{RT, String, Dict{String, Any}, Dict{String, Any}, BT, LenUnits, TD}
const ReTD{RT, BT, LenUnits} = ReT{RT, Dict{String, Any}, BT, LenUnits}
const BinaryRootedTree = ReTD{OneRoot, BinaryBranching, Float64}
const RootedTree = ReTD{OneRoot, PolytomousBranching, Float64}
const ManyRootTree = ReTD{ManyRoots, PolytomousBranching, Float64}
const UnrootedTree = ReTD{Unrooted, PolytomousBranching, Float64}

# Types matches

_emptydata(::Type{Data}) where Data = Data()
_emptydata(::Type{RecursiveNode{RT, NL, NodeData, BranchData}}) where
    {RT, NL, NodeData, BranchData} = _newdata(NodeData)
_emptydata(::Type{RecursiveBranch{RT, NL, NodeData, BranchData}}) where
    {RT, NL, NodeData, BranchData} = _newdata(BranchData)

_tiplabeltype(::Type{Nothing}) = Nothing
_tiplabeltype(::Type{<: Dict{NL}}) where NL = NL
_tiplabeltype(::Type{RecursiveTree{RT, NL, NodeData, BranchData, BT, LenUnits, TD}}) where
    {RT, NL, NodeData, BranchData, BT, LenUnits, TD} = _tiplabeltype(TD)

_matchlabels(::Type{S}, ::Type{T}) where {S, T} = false
_matchlabels(::Type{S}, ::Type{S}) where S = true
_matchlabels(::Type{S}, ::Type{Nothing}) where S = true

# Retrieving trees

import Phylo.API: _treenametype
_treenametype(::Type{<: RecursiveTree}) = String

import Phylo.API: _gettreename
_gettreename(tree::RecursiveTree) = tree.name

# Information about leaves in a single store on the tree

import Phylo.API: _leafinfotype
_leafinfotype(::Type{<:RecursiveTree{RT, NL, NodeData, BranchData, BT, LenUnits, TD}}) where
    {RT, NL, NodeData, BranchData, BT, LenUnits, TD} = TD

import Phylo.API: _getleafinfo
_getleafinfo(tree::RecursiveTree) = tree.tipdata

import Phylo.API: _setleafinfo!
_setleafinfo!(tree::RecursiveTree{RT, NL, NodeData, BranchData, BT, LenUnits, TD}, leafinfo::TD) where
    {RT, NL, NodeData, BranchData, BT, LenUnits, TD} = (tree.tipdata = leafinfo)

# Retrieving nodes

import Phylo.API: _getroots
_getroots(tree::RecursiveTree{<: Rooted}) = tree.roots

import Phylo.API: _nnodes
_nnodes(tree::RecursiveTree) = count((!)∘ismissing, tree.nodes)

import Phylo.API: _getnodes
_getnodes(tree::RecursiveTree) = skipmissing(tree.nodes)

import Phylo.API: _getnodenames
_getnodenames(tree::RecursiveTree) = keys(tree.nodedict)

import Phylo.API: _hasnode
_hasnode(tree::RecursiveTree{RT, NL}, name::NL) where {RT, NL} =
    haskey(tree.nodedict, name)
_hasnode(tree::RecursiveTree{RT, NL}, node::N) where
    {RT, NL, N <: RecursiveNode{RT, NL}} =
    haskey(tree.nodedict, node.name)

import Phylo.API: _getnode
_getnode(::RecursiveTree, node::RecursiveNode) = node
_getnode(tree::RecursiveTree{RT, NL}, name::NL) where {RT, NL} = tree.nodes[tree.nodedict[name]]

import Phylo.API: _getnodename
_getnodename(::RecursiveTree, node::RecursiveNode) = node.name
_getnodename(::RecursiveTree{RT, NL}, name::NL) where {RT, NL} = name

# Creating and destroy nodes

import Phylo.API: _createnode!
function _createnode!(tree::RecursiveTree{RT, NL, NodeData, BranchData, BT, LenUnits, TD},
                      name::Union{NL, Missing}, data::NodeData = _emptydata(NodeData)) where
    {RT <: Rooted, NL, NodeData, BranchData, BT, LenUnits, TD}

    NT = RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits}
    nodename = ismissing(name) ? _newnodelabel(tree) : name
    _hasnode(tree, nodename) && error("Node $nodename already exists in tree.")
    id = length(tree.nodes) + 1
    node = NT(name = nodename, id = id, data = data)
    push!(tree.nodes, node)
    tree.nodedict[nodename] = id
    push!(tree.roots, node)

    _invalidate!(tree)
    return node
end

function _createnode!(tree::RecursiveTree{Unrooted, NL, NodeData, BranchData, BT, LenUnits, TD},
                      name::Union{NL, Missing}, data::NodeData = _emptydata(NodeData)) where
    {NL, NodeData, BranchData, BT, LenUnits, TD}

    NT = RecursiveNode{Unrooted, NL, NodeData, BranchData, BT, LenUnits}
    nodename = ismissing(name) ? _newnodelabel(tree) : name
    _hasnode(tree, nodename) && error("Node $nodename already exists in tree.")
    id = length(tree.nodes) + 1
    node = NT(name = nodename, id = id, data = data)
    push!(tree.nodes, node)
    tree.nodedict[nodename] = id

    _invalidate!(tree)
    return node
end

import Phylo.API: _deletenode!
function _deletenode!(tree::RecursiveTree, node::RecursiveNode)
    (haskey(tree.nodedict, node.name) &&
     tree.nodedict[node.name] == node.id &&
     tree.nodes[node.id] ≡ node) ||
        error("Node $(node.name) is not in tree $(tree.name), cannot be deleted")
    
    # Does nothing for unrooted tree, as inbound is never set
    !isnothing(node.in) && _deletebranch!(tree, node.in)

    # Delete outbound connections of rooted tree, all connections of unrooted
    while _degree(tree, node) > 0
        _deletebranch!(tree, first(node.conns))
    end

    tree.nodes[tree.nodedict[node.name]] = missing
    delete!(tree.nodedict, node.name)
    filter!(n -> n ≢ node, tree.roots)
    _invalidate!(tree)
    return true
end

# Retrieving and editing connections on nodes

import Phylo.API: _hasinbound
_hasinbound(::RecursiveTree, node::RecursiveNode{<: Rooted}) = !isnothing(node.in)

import Phylo.API: _degree
_degree(::RecursiveTree, node::RecursiveNode{Unrooted}) = length(node.conns)

import Phylo.API: _getinbound
_getinbound(::RecursiveTree, node::RecursiveNode{<: Rooted}) = node.in

import Phylo.API: _getoutbounds
_getoutbounds(::RecursiveTree, node::RecursiveNode{<: Rooted}) = node.conns

import Phylo.API: _getconnections
_getconnections(::RecursiveTree, node::RecursiveNode{Unrooted}, exclude) =
    filter(∉(exclude), node.conns)

import Phylo.API: _addinbound!
function _addinbound!(tree::RecursiveTree{RT},
                      node::RecursiveNode{RT},
                      branch::RecursiveBranch{RT}) where RT <: Rooted
    _hasinbound(tree, node) &&
        error("RecursiveNode $(node.name) already has an inbound connection")
    node.in = branch
    filter!(n -> n ≠ node, tree.roots)
    _invalidate!(tree)
end

import Phylo.API: _removeinbound!
function _removeinbound!(tree::RecursiveTree{RT},
                         node::RecursiveNode{RT},
                         branch::RecursiveBranch{RT} = node.in) where RT <: Rooted
    _hasinbound(tree, node) || error("Node $(node.name) has no inbound connection")
    node.in ≡ branch ||
        error("Node $(node.name) has no inbound connection from branch $(branch.id)")
    node.in = nothing
    push!(tree.roots, node)
    _invalidate!(tree)
end

import Phylo.API: _hasoutboundspace
_hasoutboundspace(::RecursiveTree{<: Rooted, NL, TheirData, MyData,
                  BinaryBranching, LenUnits, TD},
                  node::RecursiveNode{<: Rooted, NL, TheirData, MyData,
                                      BinaryBranching, LenUnits}) where
    {NL, TheirData, MyData, TD, LenUnits} = length(node.conns) < 2

import Phylo.API: _hasspace
_hasspace(::RecursiveTree{Unrooted, NL, TheirData, MyData,
                          BinaryBranching, LenUnits, TD},
node::RecursiveNode{Unrooted, NL, TheirData, MyData,
                    BinaryBranching, LenUnits}) where
    {NL, TheirData, MyData, LenUnits, TD} = length(node.conns) < 3

import Phylo.API: _addoutbound!
function _addoutbound!(tree::RecursiveTree{RT},
                       node::RecursiveNode{RT},
                       branch::RecursiveBranch{RT}) where RT <: Rooted
    _hasoutboundspace(tree, node) ||
        error("Node $(from.name) has no outbound space")

    push!(node.conns, branch)
    _invalidate!(tree)
end

import Phylo.API: _removeoutbound!
function _removeoutbound!(tree::RecursiveTree{RT},
                          node::RecursiveNode{RT},
                          branch::RecursiveBranch{RT}) where RT <: Rooted
    if branch ∉ _getoutbounds(tree, node)
        error("Node $(node.name) does not have outbound connection to branch $(branch.id)")
    end
    filter!(p -> p ≢ branch, node.conns)
    _invalidate!(tree)
end

import Phylo.API: _addconnection!
function _addconnection!(tree::RecursiveTree{Unrooted},
                         node::RecursiveNode{Unrooted},
                         branch::RecursiveBranch{Unrooted})
    if !_hasspace(tree, node)
        error("Node $(node.name) does not have space for a new connection")
    end
    push!(node.conns, branch)
    _invalidate!(tree)
end

import Phylo.API: _removeconnection!
function _removeconnection!(tree::RecursiveTree{Unrooted},
                            node::RecursiveNode{Unrooted},
                            branch::RecursiveBranch{Unrooted})
    if branch ∉ getconnections(tree, node)
        error("Node $(node.name) does not have connection to branch $(branch.id)")
    end
    filter!(p -> p ≢ branch, node.conns)
    _invalidate!(tree)
end

# Retrieving branches

import Phylo.API: _nbranches
_nbranches(tree::RecursiveTree) = count((!)∘ismissing, tree.branches)

import Phylo.API: _getbranches
_getbranches(tree::RecursiveTree) = skipmissing(tree.branches)

import Phylo.API: _getbranchnames
_getbranchnames(tree::RecursiveTree) = [b.id for b in skipmissing(tree.branches)]

import Phylo.API: _hasbranch
_hasbranch(tree::RecursiveTree, id::Int) =
    1 ≤ id ≤ length(tree.branches) && !ismissing(tree.branches[id])
_hasbranch(tree::RecursiveTree, branch::RecursiveBranch) =
    branch ≡ tree.branches[branch.id]

import Phylo.API: _getbranch
_getbranch(tree::RecursiveTree, id::Int) = tree.branches[id]

import Phylo.API: _getbranchname
_getbranchname(::RecursiveTree, branch::RecursiveBranch) = branch.id

# Branch length info

import Phylo.API: _branchdims
_branchdims(T::Type{<:RecursiveTree}) = _branchdims(branchtype(T))
_branchdims(::Type{<:RecursiveBranch{RT, NL, NodeData, BranchData, BT, LenUnits}}) where
    {RT, NL, NodeData, BranchData, BT, LenUnits} = dimension(LenUnits)

import Phylo.API: _getlength
_getlength(::RecursiveTree, branch::RecursiveBranch) = branch.length

# Retrieving connections on branches

import Phylo.API: _src
_src(::RecursiveTree, branch::RecursiveBranch{<:Rooted}) = branch.in

import Phylo.API: _dst
_dst(::RecursiveTree, branch::RecursiveBranch{<:Rooted}) = branch.conns[1]

import Phylo.API: _conn
_conn(::RecursiveTree{Unrooted}, branch::RecursiveBranch{Unrooted},
    exclude::RecursiveNode{Unrooted}) =
        exclude ≡ branch.conns[1] ? branch.conns[2] :
            (exclude ≡ branch.conns[2] ? branch.conns[1] :
             error("Branch $(branch.id) not connected to node $(exclude.name)"))

import Phylo.API: _conns
_conns(::RecursiveTree{Unrooted}, branch::RecursiveBranch{Unrooted}) = branch.conns
         
# Information about individual nodes stored on the nodes

import Phylo.API: _nodedatatype
_nodedatatype(::Type{<:RecursiveTree{RT, NL, NodeData, BranchData}}) where
    {RT, NL, NodeData, BranchData} = NodeData
_nodedatatype(::Type{<:RecursiveNode{RT, NL, NodeData, BranchData}}) where
    {RT, NL, NodeData, BranchData} = NodeData

import Phylo.API: _getnodedata
_getnodedata(::RecursiveTree, node::RecursiveNode) = node.data

import Phylo.API: _setnodedata!
_setnodedata!(::RecursiveTree{RT, NL, NodeData, BranchData},
              node::RecursiveNode{RT, NL, NodeData, BranchData},
              data::NodeData) where {RT, NL, NodeData, BranchData} = (node.data = data)

# Information about individual branches stored on the branches

import Phylo.API: _branchdatatype
_branchdatatype(::Type{<:RecursiveTree{RT, NL, NodeData, BranchData}}) where
    {RT, NL, NodeData, BranchData} = BranchData
_branchdatatype(::Type{<:RecursiveBranch{RT, NL, NodeData, BranchData}}) where
    {RT, NL, NodeData, BranchData} = BranchData

import Phylo.API: _getbranchdata
_getbranchdata(::RecursiveTree, branch::RecursiveBranch) = branch.data

import Phylo.API: _setbranchdata!
_setbranchdata!(::RecursiveTree{RT, NL, NodeData, BranchData},
                branch::RecursiveBranch{RT, NL, NodeData, BranchData},
                data::BranchData) where {RT, NL, NodeData, BranchData} =
    (branch.data = data)

# New label methods

import Phylo.API: _newnodelabel
function _newnodelabel(tree::RecursiveTree{RT, String}) where RT
    id = length(tree.nodes) + 1
    name = NODENAME * " $id"
    return haskey(tree.nodedict, name) ?
        _newlabel(collect(getnodenames(tree)), NODENAME) : name
end

import Phylo.API: _newbranchlabel
_newbranchlabel(tree::RecursiveTree) = length(tree.branches) + 1

import Phylo.API: _createbranch!
function _createbranch!(tree::RecursiveTree{RT, NL, NodeData, BranchData, BT, LenUnits, TD},
                        from::RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits},
                        to::RecursiveNode{RT, NL, NodeData, BranchData, BT, LenUnits},
                        len::Union{LenUnits, Missing} = missing,
                        name = missing,
                        data::BranchData = _emptydata(BranchData)) where
    {RT <: Rooted, NL, NodeData, BranchData, BT, LenUnits, TD}

    if ismissing(name)
        name = length(tree.branches) + 1
    end

    ismissing(len) || len ≥ zero(len) ||
        error("Branch length must be positive or missing (no recorded length), not $len")

    _hasinboundspace(tree, to) ||
        error("Node $(to.name) already has an inbound connection")

    _hasoutboundspace(tree, from) ||
        error("Node $(from.name) has no outbound space")

    branch = RecursiveBranch{RT, NL, NodeData, BranchData,
                             BT, LenUnits}(id = name, in = from, conns = [to],
                                           length = len, data = data)
    
    if 1 ≤ name ≤ length(tree.branches)
        tree.branches[name] = branch
    else
        l = length(tree.branches) + 1
        if name > l
            resize!(tree.branches, name)
            tree.branches[l:name-1] .= missing
            tree.branches[name] = branch
        else
            push!(tree.branches, branch)
        end
    end
                                
    _addoutbound!(tree, from, branch)
    _addinbound!(tree, to, branch)

    _invalidate!(tree)
    return branch
end

function _createbranch!(tree::RecursiveTree{Unrooted, NL, NodeData, BranchData, BT, LenUnits, TD},
                        from::RecursiveNode{Unrooted, NL, NodeData, BranchData, BT, LenUnits},
                        to::RecursiveNode{Unrooted, NL, NodeData, BranchData, BT, LenUnits},
                        len::Union{LenUnits, Missing} = missing,
                        name = missing,
                        data::BranchData = _emptydata(BranchData)) where
    {Unrooted, NL, NodeData, BranchData, BT, LenUnits, TD}

    if ismissing(name)
        name = length(tree.branches) + 1
    end

    ismissing(len) || len ≥ zero(len) ||
        error("Branch length must be positive or missing (no recorded length), not $length")

    branch = RecursiveBranch{Unrooted, NL, NodeData, BranchData,
                             BT, LenUnits}(id = name, conns = [to, from],
                                           length = len, data = data)

    _hasspace(tree, from) ||
        error("Node $(from.name) has no space for new connections")

    _hasspace(tree, to) ||
        error("Node $(to.name) has no space for new connections")

    if 1 ≤ name ≤ length(tree.branches)
        tree.branches[name] = branch
    else
        l = length(tree.branches) + 1
        if name > l
            resize!(tree.branches, name)
            tree.branches[l:name-1] .= missing
            tree.branches[name] = branch
        else
            push!(tree.branches, branch)
        end
    end

    _addconnection!(tree, from, branch)
    _addconnection!(tree, to, branch)

    _invalidate!(tree)
    return branch
end

function _deletebranch!(tree::RecursiveTree{RT}, branch::RecursiveBranch{RT}) where RT <: Rooted
    _removeoutbound!(tree, branch.in, branch)
    _removeinbound!(tree, branch.conns[1], branch)
    tree.branches[branch.id] = missing
    _invalidate!(tree)
    return true
end

function _deletebranch!(tree::RecursiveTree{Unrooted}, branch::RecursiveBranch{Unrooted})
    _removeconnection!(tree, branch.conns[1], branch)
    _removeconnection!(tree, branch.conns[2], branch)
    tree.branches[branch.id] = missing
    _invalidate!(tree)
    return true
end

import Base.show
show(io::IO, ::MIME"text/plain", node::RecursiveNode{Unrooted}) =
    length(node.conns) == 0 ? print(io, "unattached node '$(node.name)'") :
        length(node.conns) == 1 ? print(io, "leaf node '$(node.name)'") :
            print(io, "internal node '$(node.name)'")

function show(io::IO, node::RecursiveNode{Unrooted})
    nc = length(node.conns)
    print(io, "RecursiveNode{Unrooted} '$(node.name)' ")
    if nc == 0
        print(io, "with no connections")
    elseif nc == 1
        print(io, "with 1 connection (branch $(node.conns[1].id))")
    else
        print(io, "with $nc outbound connections (branches $(getfield.(node.conns, :id)))")
    end
end

show(io::IO, ::MIME"text/plain", node::RecursiveNode{<: Rooted}) =
    isnothing(node.in) ?
        (length(node.conns) == 0 ? print(io, "unattached node '$(node.name)'") :
            print(io, "root node '$(node.name)'")) :
        (length(node.conns) == 0 ? print(io, "internal node '$(node.name)'") :
            print(io, "leaf node '$(node.name)'"))

function show(io::IO, node::RecursiveNode{RT}) where RT <: Rooted
    print(io, "RecursiveNode{$RT} '$(node.name)', ")
    no = length(node.conns)
    if isnothing(node.in)
        if no == 0
            print(io, "an isolated node with no connections")
        elseif no == 1
            print(io, "a root node with 1 outbound connection" *
                  " (branch $(node.conns[1].id))")
        else
            print(io, "a root node with $no outbound connections" *
                  " (branches $(getfield.(node.conns, :id)))")
        end
    else
        if no == 0
            print(io, "a leaf with an incoming connection (branch $(node.in.id))")
        elseif no == 1
            print(io, "an internal node with 1 inbound and 1 outbound connection " *
                  "(branches $(node.in.id) and $(node.conns[1].id))")
        else
            print(io, "an internal node with 1 inbound and $no outbound connections" *
                  " (branches $(node.in.id) and $(getfield.(node.conns, :id)))")
        end
    end
end

show(io::IO, ::MIME"text/plain", branch::RecursiveBranch) =
    print(io, "branch '$(branch.id)'")

function show(io::IO, branch::RecursiveBranch{Unrooted})
    print(io, "RecursiveBranch{Unrooted} $(branch.id), connecting nodes" *
          " '$(branch.conns[1].name)' and '$(branch.conns[2].name)'" *
          (ismissing(branch.length) ? "" : " (length $(branch.length))"))
end

function show(io::IO, branch::RecursiveBranch{RT}) where RT <: Rooted
    print(io, "RecursiveBranch{$RT} $(branch.id), from node '$(branch.in.name)'" *
          " to node '$(branch.conns[1].name)'" *
          (ismissing(branch.length) ? "" : " (length $(branch.length))"))
end
