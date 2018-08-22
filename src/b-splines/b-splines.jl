export
    BSpline,

    NoInterp,
    Constant,
    Linear,
    Quadratic,
    Cubic

abstract type Degree{N} <: Flag end

struct BSpline{D<:Degree} <: InterpolationType end
BSpline(::D) where {D<:Degree} = BSpline{D}()

bsplinetype(::Type{BSpline{D}}) where {D<:Degree} = D
bsplinetype(::BS) where {BS<:BSpline} = bsplinetype(BS)

degree(::BSpline{D}) where D<:Degree = D()
degree(::NoInterp) = NoInterp()

struct BSplineInterpolation{T,N,TCoefs<:AbstractArray,IT<:DimSpec{BSpline},GT<:DimSpec{GridType},Axs<:Tuple{Vararg{AbstractUnitRange,N}}} <: AbstractInterpolation{T,N,IT,GT}
    coefs::TCoefs
    parentaxes::Axs
    it::IT
    gt::GT
end
function BSplineInterpolation(::Type{TWeights}, A::AbstractArray{Tel,N}, it::IT, gt::GT, axs) where {N,Tel,TWeights<:Real,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    # String interpolation causes allocation, noinline avoids that unless they get called
    @noinline err_concrete(IT) = error("The b-spline type must be a concrete type (was $IT)")
    @noinline warn_concrete(A) = @warn("For performance reasons, consider using an array of a concrete type (typeof(A) == $(typeof(A)))")

    isconcretetype(IT) || err_concrete(IT)
    isconcretetype(typeof(A)) || warn_concrete(A)

    # Compute the output element type when positions have type TWeights
    if isempty(A)
        T = Base.promote_op(*, TWeights, eltype(A))
    else
        T = typeof(zero(TWeights) * first(A))
    end
    BSplineInterpolation{T,N,typeof(A),IT,GT,typeof(axs)}(A, fix_axis.(axs), it, gt)
end

coefficients(itp::BSplineInterpolation) = itp.coefs
interpdegree(itp::BSplineInterpolation) = interpdegree(itpflag(itp))
interpdegree(::BSpline{T}) where T = T()
interpdegree(it::Tuple{Vararg{Union{BSpline,NoInterp},N}}) where N = interpdegree.(it)
itpflag(itp::BSplineInterpolation) = itp.it
gridflag(itp::BSplineInterpolation) = itp.gt

size(itp::BSplineInterpolation) = map(length, itp.parentaxes)
axes(itp::BSplineInterpolation) = itp.parentaxes

lbounds(itp::BSplineInterpolation) = _lbounds(itp.parentaxes, itpflag(itp), gridflag(itp))
ubounds(itp::BSplineInterpolation) = _ubounds(itp.parentaxes, itpflag(itp), gridflag(itp))
_lbounds(axs, itp, gt) = (lbound(axs[1], getfirst(itp), getfirst(gt)), _lbounds(Base.tail(axs), getrest(itp), getrest(gt))...)
_ubounds(axs, itp, gt) = (ubound(axs[1], getfirst(itp), getfirst(gt)), _ubounds(Base.tail(axs), getrest(itp), getrest(gt))...)
_lbounds(::Tuple{}, itp, gt) = ()
_ubounds(::Tuple{}, itp, gt) = ()

# The unpadded defaults
lbound(ax::AbstractUnitRange, ::BSpline, ::OnCell) = first(ax) - 0.5
ubound(ax::AbstractUnitRange, ::BSpline, ::OnCell) = last(ax) + 0.5
lbound(ax::AbstractUnitRange, ::BSpline, ::OnGrid) = first(ax)
ubound(ax::AbstractUnitRange, ::BSpline, ::OnGrid) = last(ax)

fix_axis(r::Base.OneTo) = r
fix_axis(r::Base.Slice) = r
fix_axis(r::UnitRange) = Base.Slice(r)
fix_axis(r::AbstractUnitRange) = fix_axis(UnitRange(r))

count_interp_dims(::Type{BSI}, n) where BSI<:BSplineInterpolation = count_interp_dims(itptype(BSI), n)

function interpolate(::Type{TWeights}, ::Type{TC}, A, it::IT, gt::GT) where {TWeights,TC,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    Apad = prefilter(TWeights, TC, A, it, gt)
    BSplineInterpolation(TWeights, Apad, it, gt, axes(A))
end

"""
    itp = interpolate(A, interpmode, gridstyle)

Interpolate an array `A` in the mode determined by `interpmode` and `gridstyle`.
`interpmode` may be one of

- `BSpline(NoInterp())`
- `BSpline(Linear())`
- `BSpline(Quadratic(BC()))` (see [`BoundaryCondition`](@ref))
- `BSpline(Cubic(BC()))`

It may also be a tuple of such values, if you want to use different interpolation schemes along each axis.

`gridstyle` should be one of `OnGrid()` or `OnCell()`.
"""
function interpolate(A::AbstractArray, it::IT, gt::GT) where {IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    interpolate(tweight(A), tcoef(A), A, it, gt)
end

# We can't just return a tuple-of-types due to julia #12500
tweight(A::AbstractArray) = Float64
tweight(A::AbstractArray{Float32}) = Float32
tweight(A::AbstractArray{Rational{Int}}) = Rational{Int}
tweight(A::AbstractArray{T}) where {T<:Integer} = typeof(float(zero(T)))

tcoef(A::AbstractArray) = eltype(A)
tcoef(A::AbstractArray{Float32}) = Float32
tcoef(A::AbstractArray{Rational{Int}}) = Rational{Int}
tcoef(A::AbstractArray{T}) where {T<:Integer} = typeof(float(zero(T)))

function interpolate!(::Type{TWeights}, A, it::IT, gt::GT) where {TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    # Set the bounds of the interpolant inward, if necessary
    axsA = axes(A)
    axspad = padded_axes(axsA, it)
    BSplineInterpolation(TWeights, prefilter!(TWeights, A, it, gt), it, gt, fix_axis.(padinset.(axsA, axspad)))
end
function interpolate!(A::AbstractArray, it::IT, gt::GT) where {IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    interpolate!(tweight(A), A, it, gt)
end

lut!(dl, d, du) = lu!(Tridiagonal(dl, d, du), Val(false))

include("constant.jl")
include("linear.jl")
include("quadratic.jl")
include("cubic.jl")
include("indexing.jl")
include("prefiltering.jl")
include("../filter1d.jl")

Base.parent(A::BSplineInterpolation{T,N,TCoefs,UT}) where {T,N,TCoefs,UT<:Union{BSpline{Linear},BSpline{Constant}}} = A.coefs
Base.parent(A::BSplineInterpolation{T,N,TCoefs,UT}) where {T,N,TCoefs,UT} =
    throw(ArgumentError("The given BSplineInterpolation does not serve as a \"view\" for a parent array. This would only be true for Constant and Linear b-splines."))
