module Interpolations

export
    interpolate,
    interpolate!,
    extrapolate,

    gradient!,

    OnCell,
    OnGrid,

    Flat,
    Line,
    Free,
    Periodic,
    Reflect,
    Natural,
    InPlace,
    InPlaceQ

    # see the following files for further exports:
    # b-splines/b-splines.jl
    # extrapolation/extrapolation.jl

using WoodburyMatrices, Ratios, AxisAlgorithms

import Base: convert, size, getindex, gradient, promote_rule

abstract InterpolationType
immutable NoInterp <: InterpolationType end
abstract GridType
immutable OnGrid <: GridType end
immutable OnCell <: GridType end

typealias DimSpec{T} Union{T,Tuple{Vararg{Union{T,NoInterp}}},NoInterp}

abstract AbstractInterpolation{T,N,IT<:DimSpec{InterpolationType},GT<:DimSpec{GridType}} <: AbstractArray{T,N}
abstract AbstractExtrapolation{T,N,ITPT,IT,GT} <: AbstractInterpolation{T,N,IT,GT}

abstract BoundaryCondition
immutable Flat <: BoundaryCondition end
immutable Line <: BoundaryCondition end
immutable Free <: BoundaryCondition end
immutable Periodic <: BoundaryCondition end
immutable Reflect <: BoundaryCondition end
immutable InPlace <: BoundaryCondition end
# InPlaceQ is exact for an underlying quadratic. This is nice for ground-truth testing
# of in-place (unpadded) interpolation.
immutable InPlaceQ <: BoundaryCondition end
typealias Natural Line

# TODO: size might have to be faster?
size{T,N}(itp::AbstractInterpolation{T,N}) = ntuple(i->size(itp,i), N)::NTuple{N,Int}
size(exp::AbstractExtrapolation, d) = size(exp.itp, d)
 itptype{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}) = IT
gridtype{T,N,IT,GT}(itp::AbstractInterpolation{T,N,IT,GT}) = GT

@inline gradient{T,N}(itp::AbstractInterpolation{T,N}, xs...) = gradient!(Array(T,N), itp, xs...)

include("nointerp/nointerp.jl")
include("b-splines/b-splines.jl")
include("gridded/gridded.jl")
include("extrapolation/extrapolation.jl")

end # module
