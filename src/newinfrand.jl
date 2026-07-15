
############################################################
# BrownianTipAnyDist
############################################################

using Random
using Distributions
using ReverseDiff
using Bijectors
using Phylo
using LinearAlgebra

import Distributions: rand, logpdf
import Bijectors

# ============================================================
# Helpers
# ============================================================

_untrack(x) = x
_untrack(x::ReverseDiff.TrackedReal) = ReverseDiff.value(x)
_untrack(x::ReverseDiff.TrackedArray) = ReverseDiff.value(x)

@inline _valid_λ(λ) = isfinite(λ) && (0 < λ < 1)

@inline function _valid_cov(Σ)
    any(!isfinite, Σ) && return false
    F = cholesky(Hermitian(Σ); check=false)
    return all(>(0), real.(diag(F.U)))
end

function loglik2(n, nd, σ², β)
    return -(1/2) * (
        n * log(2π) +
        nd.logV +
        (1/σ²) * (nd.yy[] - 2*nd.Q[]*β + nd.xx * β^2)
    )
end

# ============================================================
# SPD Cholesky bijector
# ============================================================

struct SPDCholeskyBijector
    m::Int
end

@inline _npacked(m::Int) = m*(m+1) ÷ 2

function Bijectors.transform(b::SPDCholeskyBijector, η::AbstractVector)
    m = b.m
    L = zeros(eltype(η), m, m)
    k = 1
    @inbounds for j in 1:m
        for i in j:m
            L[i,j] = (i==j) ? exp(η[k]) : η[k]
            k += 1
        end
    end
    return L*L'
end

function Bijectors.inv(b::SPDCholeskyBijector, Σ::AbstractMatrix)
    m = b.m
    F = cholesky(Hermitian(Σ); check=true)
    L = F.L
    η = similar(vec(L), _npacked(m))
    k = 1
    @inbounds for j in 1:m
        for i in j:m
            η[k] = (i==j) ? log(L[i,j]) : L[i,j]
            k += 1
        end
    end
    return η
end

function Bijectors.logabsdetjac(b::SPDCholeskyBijector, η::AbstractVector)
    m = b.m
    logJ = m*log(2)
    k = 1
    for j in 1:m
        logJ += (m+2-j)*η[k]
        k += (m-j+1)
    end
    return logJ
end

λ_bijector() = Bijectors.Logit()

# ============================================================
# Cache (structural only)
# ============================================================

mutable struct BrownianTipAnyCache{T<:AbstractTree}
    master_tree::T
    nodes_post::Vector
    leaves::Vector
    m::Int
    trait::Vector{String}

    leaf_len::Vector{Float64}
    leaf_ybuf::Vector{Vector{Any}}  
    len_by_node::Vector{Float64}
    htor_by_node::Vector{Float64}
end

function BrownianTipAnyCache(tree::T, traitnames::Vector{String}) where {T<:AbstractTree}

    nodeT = nodetype(typeof(tree))

    nodes_post = Vector{nodeT}(getnodes(tree, postorder))
    leaves     = Vector{nodeT}(getleaves(tree, postorder))
    m          = length(traitnames)

    leaf_len  = Float64[getlength(tree, getinbound(tree, leaf)) for leaf in leaves]

    # Vector{Any} but initially populated with zeros to avoid undef
    leaf_ybuf = [Any[0.0 for _ in 1:m] for _ in leaves]

    len_by_node = Float64[
        isroot(tree,node) ? 0.0 : getlength(tree, getinbound(tree,node))
        for node in nodes_post
    ]

    htor_by_node = Float64[
        isroot(tree,node) ? 0.0 : heighttoroot(tree,node)
        for node in nodes_post
    ]

    # initial dummy traitdata (required by Phylo)
    for (i,leaf) in enumerate(leaves)
        td = traitdata(eltype(nodedatatype(typeof(tree))),
                       traitnames,
                       leaf_ybuf[i],
                       leaf_len[i])
        setnodedata!(tree, leaf, td)
    end

    return BrownianTipAnyCache(
        tree, nodes_post, leaves, m, traitnames,
        leaf_len, leaf_ybuf, len_by_node, htor_by_node
    )
end

# ============================================================
# Distribution Type
# ============================================================

mutable struct BrownianTipAnyDist{T<:AbstractTree, N<:Number} <: ContinuousMultivariateDistribution
    Σ::Matrix{N}
    β::AbstractVector{N}
    λ::Union{Nothing, N, ReverseDiff.TrackedReal}
    tree::T
    f::Function
    cache::BrownianTipAnyCache{T}
end


Base.length(d::BrownianTipAnyDist) = nleaves(d.tree) * d.cache.m
Base.size(d::BrownianTipAnyDist) = (length(d),)
Base.size(d::BrownianTipAnyDist, i::Integer) = (i == 1 ? length(d) : 1)

# ============================================================
# Constructors
# ============================================================

function BrownianTipAnyDist(tree::T; β, Σ, trait=nothing, λ=nothing,
                            f::Function=identity) where {T<:AbstractTree}

    traitnames = trait === nothing ?
        getnodedata(tree, getroot(tree)).name :
        trait

    traitnames = isa(traitnames,String) ? [traitnames] : String.(traitnames)

    cache = BrownianTipAnyCache(tree, traitnames)
    return BrownianTipAnyDist{T,eltype(Σ)}(Σ, β, λ, tree, f, cache)
end

function BrownianTipAnyDist(tree::T, trait::AbstractString; β, σ², λ=nothing,
                            f::Function=identity) where {T<:AbstractTree}
    return BrownianTipAnyDist(tree; β=[β], Σ=reshape([σ²],1,1),
                              trait=[trait], λ=λ, f=f)
end

# ============================================================
# Bijectors
# ============================================================

parameter_bijectors(d::BrownianTipAnyDist) =
    d.λ === nothing ?
        (Σ = SPDCholeskyBijector(d.cache.m), β=identity) :
        (Σ = SPDCholeskyBijector(d.cache.m), β=identity, λ=λ_bijector())

Bijectors.bijector(::BrownianTipAnyDist) = identity

# ============================================================
# rand 
# ============================================================

function _rand_univariate(rng, d::BrownianTipAnyDist, β, λ, σ²)
    nodeT = nodetype(typeof(d.tree))
    nameT = nodenametype(typeof(d.tree))

    untrait = Dict{nodeT, Vector{Float64}}()
    traitbyname = Dict{nameT, Any}()

    for node in traversal(d.tree, preorder)
        if isroot(d.tree,node)
            v = [β[1]]
            untrait[node] = v
            traitbyname[getnodename(d.tree,node)] = d.f(v)
            continue
        end

        inb  = getinbound(d.tree,node)
        prt  = src(d.tree,inb)
        prev = untrait[prt][1]

        blen = if λ === nothing
            getlength(d.tree,inb)
        elseif isleaf(d.tree,node)
            λ*getlength(d.tree,inb) + (1-λ)*getheight(d.tree,node)
        else
            λ*getlength(d.tree,inb)
        end

        x = prev + randn(rng) * sqrt(σ² * blen)
        v = [x]

        untrait[node] = v
        traitbyname[getnodename(d.tree,node)] = d.f(v)
    end

    leafnames = getleafnames(d.tree, postorder)
    z = Vector{Float64}(undef, length(leafnames))

    @inbounds for i in eachindex(leafnames)
        z[i] = traitbyname[leafnames[i]][1]
    end

    return z
end

function _rand_multivariate(rng, d::BrownianTipAnyDist, β, λ, L)
    nodeT = nodetype(typeof(d.tree))
    nameT = nodenametype(typeof(d.tree))

    untrait = Dict{nodeT, Vector{Float64}}()
    traitbyname = Dict{nameT, Any}()

    m = d.cache.m

    for node in traversal(d.tree, preorder)
        if isroot(d.tree,node)
            v = copy(β)
            untrait[node] = v
            traitbyname[getnodename(d.tree,node)] = d.f(v)
            continue
        end

        inb = getinbound(d.tree,node)
        prt = src(d.tree,inb)
        prev = untrait[prt]

        blen = if λ === nothing
            getlength(d.tree,inb)
        elseif isleaf(d.tree,node)
            λ*getlength(d.tree,inb) + (1-λ)*getheight(d.tree,node)
        else
            λ*getlength(d.tree,inb)
        end

        inc = L * randn(rng, m) * sqrt(blen)
        v = prev .+ inc

        untrait[node] = v
        traitbyname[getnodename(d.tree,node)] = d.f(v)
    end

    leafnames = getleafnames(d.tree, postorder)
    z = Vector{Float64}(undef, m * length(leafnames))

    k = 1
    for leaf in leafnames
        v = traitbyname[leaf]
        @inbounds for j in 1:m
            z[k] = v[j]
            k += 1
        end
    end

    return z
end

function Distributions.rand(rng::AbstractRNG, d::BrownianTipAnyDist)
    β = Vector{Float64}(_untrack(d.β))
    λ = d.λ === nothing ? nothing : Float64(_untrack(d.λ))
    m = d.cache.m

    if m == 1
        σ² = Float64(_untrack(d.Σ[1,1]))
        return _rand_univariate(rng, d, β, λ, σ²)
    end

    Σ = Matrix{Float64}(_untrack(d.Σ))
    L = cholesky(Hermitian(Σ)).L
    return _rand_multivariate(rng, d, β, λ, L)
end

# ============================================================
# logpdf 
# ============================================================

function logpdf(d::BrownianTipAnyDist, z::AbstractVector{<:Number})

    c = d.cache
    m = c.m
    n = length(c.leaves)

    tree = c.master_tree
    
    root = getroot(tree)
    getnodedata(tree, root).t = eps()   # tiny branch length

    length(z) == n*m ||
        throw(DimensionMismatch("Expected length $(n*m)"))

    Σu = _untrack(d.Σ)
    _valid_cov(Σu) || return -Inf
    if d.λ !== nothing && !_valid_λ(d.λ)
        return -Inf
    end

    @inbounds for i in 1:n
        src = (i-1)*m + 1
        for j in 1:m
            c.leaf_ybuf[i][j] = z[src + j - 1]
        end
    end

    nodes = c.nodes_post
    trait = c.trait


    # -------------------------------------------------------
    # apply λ
    # -------------------------------------------------------

    if d.λ === nothing
        @inbounds for k in eachindex(nodes)
            node = nodes[k]
            ndn  = getnodedata(tree,node)
            if isroot(tree,node)
                ndn.t = 0.0
            else
                ndn.t = c.len_by_node[k]
            end
        end
    else
        λ = d.λ
        @inbounds for k in eachindex(nodes)
            node = nodes[k]
            ndn  = getnodedata(tree,node)
            if isroot(tree,node)
                ndn.t = zero(λ)
            elseif isleaf(tree,node)
                ndn.t =
                    λ * c.len_by_node[k] +
                    (1 - λ) * c.htor_by_node[k]
            else
                ndn.t =
                    λ * c.len_by_node[k]
            end
        end
    end



    # -------------------------------------------------------
    # Detect whether z is latent or observed
    # -------------------------------------------------------
    latent_z = any(x -> x isa ReverseDiff.TrackedReal, z)

    # Always overwrite traitdata for BM, univariate or multivariate
    @inbounds for i in 1:n
        td = traitdata(
            eltype(nodedatatype(typeof(tree))),
            trait,
            c.leaf_ybuf[i],
            c.leaf_len[i]
        )
        setnodedata!(tree, c.leaves[i], td)
    end


    
    # -------------------------------------------------------
    # threepoint!
    # -------------------------------------------------------
    threepoint!(tree, trait, nodes)
    nd = getnodedata(tree, last(nodes))



    A = nd.yy .- 2*(nd.Q * d.β') .+ (d.β * nd.xx * d.β')
    logdetΣ = logabsdet(Σu)[1]
    trterm  = tr(Σu \ A)

    return -(1/2)*(
        n*m*log(2π) +
        m*nd.logV +
        n*logdetΣ +
        trterm
    )
end


loglikelihood(d::BrownianTipAnyDist, z::AbstractVector) = logpdf(d,z)

