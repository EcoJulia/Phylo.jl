import Phylo.API: _hasheight, _getheight, _setheight!

mutable struct LeafInfo <: AbstractInfo
    height::Float64
end

LeafInfo() = LeafInfo(NaN)

_hasheight(li::LeafInfo) = !isnan(li.height)
_getheight(li::LeafInfo) = li.height
_setheight!(li::LeafInfo, height::Float64) = (li.height = height)

#  - _hasheight()
function _hasheight(::AbstractInfo)
    return false
end

#  - _getheight() - should never be called as hasheight returns false
function _getheight(::AbstractInfo)
    throw(NullException())
    return NaN
end

#  - _setheight!() - ignore set value
function _setheight!(::AbstractInfo, value::Float64)
    warn("Ignoring height")
    return value
end
