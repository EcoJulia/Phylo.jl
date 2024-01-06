using RecipesBase

@recipe function f(tree::Phylo.AbstractTree; treetype = :dendrogram,
                   marker_group = nothing, line_group = nothing,
                   showtips = true, tipfont = (7,))
    linecolor --> :black
    grid --> false
    framestyle --> :none
    legend --> false
    colorbar --> true
    size --> (1000, 1000)

    d, h, n = _findxy(tree)
    adj = 0.03maximum(values(d))
    tipannotations = map(x -> (d[x] + adj, h[x], x), getleafnames(tree))

    lz = get(plotattributes, :line_z, nothing)
    mz = get(plotattributes, :marker_z, nothing)
    isnothing(lz) || (line_z := _handlez(lz, tree, n))
    isnothing(mz) || (marker_z := _handlez(mz, tree, n))
    mg = _handlez(marker_group, tree, n)
    lg = _handlez(line_group, tree, n)

    x, y = Float64[], Float64[]
    for node in n
        if hasinbound(tree, node)
            push!(x, d[getparent(tree, node)], d[getparent(tree, node)],
                  d[node], NaN)
            push!(y, h[getparent(tree, node)], h[node], h[node], NaN)
        end
    end

    marker_x, marker_y = _handlemarkers(plotattributes, mg, tree, d, h, n)

    if treetype == :dendrogram
        Dendrogram(x, y, tipannotations, marker_x, marker_y, showtips, tipfont,
                   mg, lg)
    elseif treetype == :fan
        Fan(x, y, tipannotations, marker_x, marker_y, showtips, tipfont, mg, lg)
    else
        throw(ArgumentError("Unsupported `treetype`; valid values are `:dendrogram` or `:fan`"))
    end
end

struct Dendrogram
    x::Any
    y::Any
    tipannotations::Any
    marker_x::Any
    marker_y::Any
    showtips::Any
    tipfont::Any
    marker_group::Any
    line_group::Any
end
struct Fan
    x::Any
    y::Any
    tipannotations::Any
    marker_x::Any
    marker_y::Any
    showtips::Any
    tipfont::Any
    marker_group::Any
    line_group::Any
end

@recipe function f(dend::Dendrogram)
    ex = extrema(filter(isfinite, dend.x))
    xlims --> (ex[1] - 0.05 * ex[2], ex[2] * 1.15)

    # tip annotations
    dend.showtips &&
        (annotations := map(x -> (x[1], x[2], (x[3], :left, dend.tipfont...)),
                            dend.tipannotations))

    sa = get(plotattributes, :series_annotations, nothing)
    @series begin
        seriestype := :path
        markersize := 0
        markershape := :none
        series_annotations := nothing
        label --> ""
        #    primary := false

        lc = _extend(get(plotattributes, :linecolor, nothing), dend.x)
        lc !== nothing && (linecolor := lc)
        la = _extend(get(plotattributes, :linealpha, nothing), dend.x)
        la !== nothing && (linealpha := la)
        lz = _extend(get(plotattributes, :line_z, nothing), dend.x)
        lz !== nothing && (line_z := lz)

        dend.x, dend.y
    end
    if !isempty(dend.marker_x) || sa !== nothing
        if isnothing(dend.marker_group)
            @series begin
                seriestype := :scatter
                sa !== nothing && (series_annotations := sa)
                label --> ""
                dend.marker_x, dend.marker_y
            end
        else
            groups = sort(unique(dend.marker_group))
            for group in groups
                idxs = findall(==(group), dend.marker_group)
                @series begin
                    seriestype := :scatter
                    sa !== nothing && (series_annotations := sa[idxs])
                    label --> string(group)
                    dend.marker_x[idxs], dend.marker_y[idxs]
                end
            end
        end
    end
    primary := false
    label := ""
    return nothing
end

@recipe function f(fan::Fan)
    adjust(y) = 2pi * y / (length(fan.tipannotations) + 1)

    aspect_ratio := 1

    # tip annotations
    mx = maximum(filter(isfinite, fan.x))
    if fan.showtips
        xlims --> (1.5 .* (-mx, mx))
        ylims --> (1.5 .* (-mx, mx))
        annotations := map(x -> (_tocirc(x[1], adjust(x[2]))...,
                                 (x[3], :left,
                                  rad2deg(adjust(x[2])), fan.tipfont...)),
                           fan.tipannotations)
    end

    sa = get(plotattributes, :series_annotations, nothing)
    @series begin
        seriestype := :path
        markersize := 0
        markershape := :none
        series_annotations := nothing
        label := ""

        x, y = _circle_transform_segments(fan.x, adjust(fan.y))
        lc = _extend(get(plotattributes, :linecolor, nothing), x)
        lc !== nothing && (linecolor := lc)
        la = _extend(get(plotattributes, :linealpha, nothing), x)
        la !== nothing && (linealpha := la)
        lz = _extend(get(plotattributes, :line_z, nothing), x)
        lz !== nothing && (line_z := lz)
        x, y
    end
    if !isempty(fan.marker_x) || sa !== nothing
        if isnothing(fan.marker_group)
            @series begin
                seriestype := :scatter
                sa !== nothing && (series_annotations := sa)
                label --> ""
                _xcirc.(adjust(fan.marker_y), fan.marker_x),
                _ycirc.(adjust(fan.marker_y), fan.marker_x)
            end
        else
            groups = sort(unique(fan.marker_group))
            for group in groups
                idxs = findall(==(group), fan.marker_group)
                @series begin
                    seriestype := :scatter
                    sa !== nothing && (series_annotations := sa[idxs])
                    label --> string(group)
                    _xcirc.(adjust(fan.marker_y[idxs]), fan.marker_x[idxs]),
                    _ycirc.(adjust(fan.marker_y[idxs]), fan.marker_x[idxs])
                end
            end
        end
    end
    primary := false
    label := ""
    return nothing
end

_handlez(x, tree, names) = x
function _handlez(x::Union{String, Symbol}, tree, names)
    return [getnodedata(tree, name, x) for name in names]
end
function _handlez(x::Dict, tree, names)
    ret = fill(NaN, length(names))
    for (i, n) in enumerate(names)
        name = getnodename(tree, n)
        if haskey(x, name)
            ret[i] = x[name]
        end
    end
    return ret
end

_mylength(x) = 1
_mylength(x::AbstractVector) = length(x)
function _handlemarkers(plotattributes, marker_group, tree, d, h, names)
    isnothing(marker_group) || (plotattributes[:marker_group] = marker_group)
    markerfields = filter(x -> occursin(r"marker", String(x)),
                          keys(plotattributes))
    isempty(markerfields) && return (Float64[], Float64[])
    maxlength = maximum([_mylength(plotattributes[k]) for k in markerfields])
    if maxlength ∈ (1, ninternal(tree))
        names = filter(node -> !isleaf(tree, node), names)
    end
    return [d[node] for node in names], [h[node] for node in names]
end

function _extend(tmp, x)
    tmp isa AbstractVector && abs(length(tmp) - count(isnan, x)) < 2 ||
        return nothing
    ret = similar(x, eltype(tmp))
    j = 1 + length(tmp) - count(isnan, x)
    for i in eachindex(x)
        ret[i] = tmp[j]
        isnan(x[i]) && (j += 1)
    end
    return ret
end

"""
    sort!(::AbstractTree; rev = false)

Sorts the branches descending from each node by total number of
descendants. This creates a clearer tree for plotting. The
process is also called "ladderizing" the tree. Use `rev=true` to
reverse the sorting order.
"""
function Base.sort!(tree::T; rev = false) where {T <: AbstractTree}
    function loc!(clade::String)
        if isleaf(tree, clade)
            return 1
        end

        sizes = map(loc!, getchildren(tree, clade))
        node = getnode(tree, clade)
        if T <: LinkTree
            node.other .= node.other[sortperm(sizes, rev = rev)]
        elseif T <: RecursiveTree
            node.conns .= node.conns[sortperm(sizes, rev = rev)]
        end
        return sum(sizes) + 1
    end

    loc!(first(nodenamefilter(isroot, tree)))
    return tree
end

"""
    sort(::AbstractTree; rev = false)

Copies a tree and sorts its branches. See `sort!` for further details.
"""
Base.sort(tree::AbstractTree; rev = false) = sort!(copy(tree))

function _nodedepths(tree::Phylo.AbstractTree)
    function finddepths!(clade::String)
        if !in(clade, keys(depth))
            for subclade in getchildren(tree, clade)
                finddepths!(subclade)
            end
        end
        if !isleaf(tree, clade)
            ch_depths = [depth[child] for child in getchildren(tree, clade)]
            depth[clade] = (maximum(ch_depths) + minimum(ch_depths)) / 2.0
        end
    end

    depth = Dict(tip => float(i)
                 for (i, tip) in enumerate(getnodename(tree, x)
                                           for x in traversal(tree, preorder)
                                           if isleaf(tree, x)))
    sizehint!(depth, nnodes(tree))
    root = getnodename(tree, getroot(tree))
    finddepths!(root)
    return depth
end

function _findxy(tree::Phylo.AbstractTree)
    height = nodeheights(tree)
    depth = _nodedepths(tree)
    return height, depth, AxisArrays.axes(height, 1).val
end

function _find_tips(height, depth, tree)
    x, y, l = Float64[], Float64[], String[]
    for k in keys(height)
        if isleaf(tree, k)
            push!(x, height[k])
            push!(y, depth[k])
            push!(l, k)
        end
    end
    return x, y, l
end

function _p_circ(start_θ, end_θ, r = 1)
    steps = range(start_θ, stop = end_θ,
                  length = 1 + ceil(Int, 60abs(end_θ - start_θ)))
    retx = Array{Float64}(undef, length(steps))
    rety = similar(retx)
    for u in eachindex(steps)
        retx[u] = _xcirc(steps[u], r)
        rety[u] = _ycirc(steps[u], r)
    end
    return retx, rety
end

_xcirc(x, r) = r * cos(x)
_ycirc(y, r) = r * sin(y)
_tocirc(x, y) = _xcirc(y, x), _ycirc(y, x)

function _circle_transform_segments(xs, ys)
    retx, rety = Float64[], Float64[]
    function _transform_seg(_x, _y)
        tmpx, tmpy = _p_circ(_y[1], _y[2], _x[1])
        append!(retx, tmpx)
        append!(rety, tmpy)
        push!(retx, _xcirc(_y[3], _x[3]), NaN)
        return push!(rety, _ycirc(_y[3], _x[3]), NaN)
    end
    i = 1
    while !isnothing(i) && i < length(xs)
        j = findnext(isnan, xs, i) - 1
        _transform_seg(view(xs, i:j), view(ys, i:j))
        i = j + 2
    end
    return retx, rety
end

"""
    map_depthfirst(FUN, start, tree, eltype = nothing)

Apply `FUN` to each node in `tree` in depth-first order, and return the result.
`FUN` must take two arguments, `val` and `node`,  where `val` is the result of
applying `FUN` to the previous node, and `node` is the current node. `start`
specifies the initial value of `val`, and `eltype` specifies the type of the
return value of `FUN`.

### Examples

≡≡≡≡≡≡≡≡≡≡≡
Define a function to evolve a trait on the tree according to Brownian motion

julia> evolve(tree) = map_depthfirst((val, node) -> val + randn(), 0., tree, Float64)
"""
function map_depthfirst(FUN, start, tree, eltype = nothing)
    root = first(nodenamefilter(isroot, tree))
    isnothing(eltype) && (eltype = typeof(FUN(start, root)))
    ret = Vector{eltype}()
    function local!(val, node)
        push!(ret, val)
        for ch in getchildren(tree, node)
            local!(FUN(val, node), ch)
        end
    end
    local!(start, root)
    return ret
end
