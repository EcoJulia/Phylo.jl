importall Phylo.API

type LeafInfo <: AbstractInfo
    height::Nullable{Float64}
end

LeafInfo() = LeafInfo(Nullable{Float64}())

_hasheight(li::LeafInfo) = !isnull(li.height)
_getheight(li::LeafInfo) = get(li.height)
_setheight!(li::LeafInfo, height::Float64) = (li.height = height)

#  - _hasheight()
function _hasheight(::AbstractInfo)
    return false
end
function _hasheight(::Void)
    return false
end
#  - _getheight() - should never be called as hasheight returns false
function _getheight(::AbstractInfo)
    throw(NullException())
    return NaN
end
function _getheight(::Void)
    throw(NullException())
    return NaN
end
#  - _setheight!() - ignore set value
function _setheight!(::AbstractInfo, value::Float64)
    warn("Ignoring height")
    return value
end
function _setheight!(::Void, value::Float64)
    warn("Ignoring height")
    return value
end
