module Interpolations

export
    interpolate,

    OnCell,
    OnGrid,

    getindex

    # see the following files for further exports:
    # b-splines/b-splines.jl

import Base: size, getindex

abstract InterpolationType
abstract GridType
immutable OnGrid <: GridType end
immutable OnCell <: GridType end

abstract AbstractInterpolation{T,N,IT<:InterpolationType,GT<:GridType} <: AbstractArray{T,N}
abstract AbstractExtrapolation{T,N,IT,GT} <: AbstractInterpolation{T,N,IT,GT}

abstract BoundaryCondition
immutable None <: BoundaryCondition end
immutable Flat <: BoundaryCondition end
immutable Line <: BoundaryCondition end
immutable Free <: BoundaryCondition end
immutable Periodic <: BoundaryCondition end
immutable Reflecting <: BoundaryCondition end
typealias Natural Line

size(itp::AbstractInterpolation) = tuple([size(itp,i) for i in 1:ndims(itp)]...)
# Fix display issues by defining odd combinations of dimensionality and index
function getindex{T}(itp::AbstractInterpolation{T,1}, x::Real, y::Real)
    if y != 1
        throw(ArgumentError("Cannot index into 1D interpolation with two indices unless y == 1 (y == $y)"))
    end
    itp[x]
end


include("b-splines/b-splines.jl")
include("extrapolation/extrapolation.jl")

end
