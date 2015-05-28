module Interpolations

export
    interpolate,
    extrapolate,

    OnCell,
    OnGrid,

    Flat,
    Line,
    Free,
    Periodic,
    Reflect,
    Natural

    # see the following files for further exports:
    # b-splines/b-splines.jl
    # extrapolation/extrapolation.jl

using WoodburyMatrices

import Base: size, getindex

abstract InterpolationType
abstract GridType
immutable OnGrid <: GridType end
immutable OnCell <: GridType end

abstract AbstractInterpolation{T,N,IT<:InterpolationType,GT<:GridType} <: AbstractArray{T,N}
abstract AbstractExtrapolation{T,N,ITPT,IT,GT} <: AbstractInterpolation{T,N,IT,GT}

abstract BoundaryCondition
immutable Flat <: BoundaryCondition end
immutable Line <: BoundaryCondition end
immutable Free <: BoundaryCondition end
immutable Periodic <: BoundaryCondition end
immutable Reflect <: BoundaryCondition end
typealias Natural Line

size(itp::AbstractInterpolation) = tuple([size(itp,i) for i in 1:ndims(itp)]...)
size(exp::AbstractExtrapolation, d) = size(exp.itp, d)
gridtype{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}) = GT

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
