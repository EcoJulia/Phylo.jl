using RecipesBase
using Phylo
using Plots # now necessary for the labels :-(

@recipe function f(tree::Phylo.AbstractTree; treetype = :dendrogram, showtips = true, tipfont = (5,))

    linecolor --> :black
    grid --> false
    framestyle --> :none
    legend --> false
    colorbar --> true
    size --> (1000, 1000)

    d, h, n = _findxy(tree)
    adj = 0.03maximum(values(d))
    tipannotations = map(x->(d[x] + adj, h[x], x), leafiter(tree))

    x, y = Float64[], Float64[]
    for node ∈ n
        if hasinbound(tree, node)
            push!(x, d[getparent(tree, node)], d[getparent(tree, node)], d[node], NaN)
            push!(y, h[getparent(tree, node)], h[node], h[node], NaN)
        end
    end

    marker_x, marker_y = Float64[], Float64[]
    if any(x->occursin(r"marker", String(x)), keys(plotattributes))
        f = nodenamefilter(!isleaf, tree)
        append!(marker_x, getindex.(Ref(d), f))
        append!(marker_y, getindex.(Ref(h), f))
    end

    if treetype == :dendrogram
        Dendrogram(x, y, tipannotations, marker_x, marker_y, showtips, tipfont)
    elseif treetype == :fan
        Fan(x, y, tipannotations, marker_x, marker_y, showtips, tipfont)
    else
        throw(ArgumentError("Unsupported `treetype`; valid values are `:dendrogram` or `:fan`"))
    end
end

struct Dendrogram; x; y; tipannotations; marker_x; marker_y; showtips; tipfont; end
struct Fan; x; y; tipannotations; marker_x; marker_y; showtips; tipfont; end

@recipe function f(d::Dendrogram)
    @series begin
        seriestype := :path
        markersize := 0
        markershape := :none

        lc = _extend(get(plotattributes, :linecolor, nothing), d.x)
        lc !== nothing && (linecolor := lc)
        la = _extend(get(plotattributes, :linealpha, nothing), d.x)
        la !== nothing && (linealpha := la)
        lz = _extend(get(plotattributes, :line_z, nothing), d.x)
        lz !== nothing && (line_z := lz)

        d.x, d.y
    end
    if !isempty(d.marker_x)
        @series begin
            seriestype := :scatter
            d.marker_x, d.marker_y
        end
    end
    d.showtips && (annotations := map(x -> (x[1], x[2], text(x[3], :left, d.tipfont...)), d.tipannotations))
    [],[]
end

@recipe function f(d::Fan)
    adjust(y) = 2pi*y / (length(d.tipannotations) + 1)
    @series begin
        seriestype := :path
        markersize := 0
        markershape := :none
        x, y = _circle_transform_segments(d.x, adjust(d.y))
        lc = _extend(get(plotattributes, :linecolor, nothing), x)
        lc !== nothing && (linecolor := lc)
        la = _extend(get(plotattributes, :linealpha, nothing), x)
        la !== nothing && (linealpha := la)
        lz = _extend(get(plotattributes, :line_z, nothing), x)
        lz !== nothing && (line_z := lz)
        x, y
    end
    if !isempty(d.marker_x)
        @series begin
            seriestype := :scatter
            _xcirc.(adjust(d.marker_y), d.marker_x), _ycirc.(adjust(d.marker_y), d.marker_x)
        end
    end
    aspect_ratio := 1
    mx = maximum(filter(isfinite, d.x))
    if d.showtips
        xlim --> (1.3 .* (-mx, mx))
        ylim --> (1.3 .* (-mx, mx))
        annotations := map(x -> (_tocirc(x[1], adjust(x[2]))..., text(x[3], :left,
            rad2deg(adjust(x[2])), d.tipfont...)), d.tipannotations)
    end
    [],[]
end

function _extend(tmp, x)
    tmp isa AbstractVector && abs(length(tmp) - count(isnan, x)) < 2 || return nothing
    ret = similar(x, eltype(tmp))
    j = 1 + length(tmp) - count(isnan, x)
    for i in eachindex(x)
        ret[i] = tmp[j]
        isnan(x[i]) && (j += 1)
    end
    return ret
end

leafiter(tree) = nodenamefilter(isleaf, tree)

function _findxy(tree::Phylo.AbstractTree)

    # two convenience recursive functions using captured variables
    function findheights!(clade::String)
        if !in(clade, keys(height))
            for subclade in getchildren(tree, clade)
                findheights!(subclade)
            end
        end
        if !isleaf(tree, clade)
            ch_heights = [height[child] for child in getchildren(tree, clade)]
            height[clade] = (maximum(ch_heights) + minimum(ch_heights)) / 2.
        end
    end

    function finddepths!(clade::String)
        push!(names, clade)
        if hasinbound(tree, clade)
            depth[clade] = depth[getparent(tree, clade)] + getbranch(tree, getinbound(tree, clade)).length
        end
        for ch in getchildren(tree, clade)
            finddepths!(ch)
        end
    end

    root = first(nodenamefilter(isroot, tree))
    height = Dict(tip => float(i) for (i, tip) in enumerate(leafiter(tree)))
    sizehint!(height, length(nodeiter(tree)))
    findheights!(root)

    depth = Dict{String, Float64}(root => 0)
    names = String[]
    sizehint!(depth, length(nodeiter(tree)))
    sizehint!(names, length(nodeiter(tree)))
    finddepths!(root)

    depth, height, names
end

function _find_tips(depth, height, tree)
    x, y, l = Float64[], Float64[], String[]
    for k in keys(depth)
        if isleaf(tree, k)
            push!(x, depth[k])
            push!(y, height[k])
            push!(l, k)
        end
    end
    x, y, l
end

function _p_circ(start_θ, end_θ, r=1)
    steps = range(start_θ, stop=end_θ, length = 1+ceil(Int, 60abs(end_θ - start_θ)))
    retx = Array{Float64}(undef, length(steps))
    rety = similar(retx)
    for u in eachindex(steps)
        retx[u] = _xcirc(steps[u], r)
        rety[u] = _ycirc(steps[u], r)
    end
    retx, rety
end

_xcirc(x, r) = r*cos(x)
_ycirc(y, r) = r*sin(y)
_tocirc(x, y) = _xcirc(y, x), _ycirc(y, x)

function _circle_transform_segments(xs, ys)
    retx, rety = Float64[], Float64[]
    function _transform_seg(_x, _y)
        tmpx, tmpy = _p_circ(_y[1], _y[2], _x[1])
        append!(retx, tmpx)
        append!(rety, tmpy)
        push!(retx, _xcirc(_y[3], _x[3]), NaN)
        push!(rety, _ycirc(_y[3], _x[3]), NaN)
    end
    i = 1
    while !(i === nothing) && i < length(xs)
        j = findnext(isnan, xs, i) - 1
        _transform_seg(view(xs,i:j), view(ys, i:j))
        i = j + 2
    end
    retx, rety
end
