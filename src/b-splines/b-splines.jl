export
    BSpline,

    NoInterp,
    Constant,
    Linear,
    Quadratic,
    Cubic

abstract type Degree{N} <: Flag end
abstract type DegreeBC{N} <: Degree{N} end  # degree type supporting a BoundaryCondition

"""
    BSpline(degree)

A flag signaling `BSpline` (integer-grid b-spline) interpolation along the corresponding axis.
`degree` is one of [`Constant`](@ref), [`Linear`](@ref), [`Quadratic`](@ref), or [`Cubic`](@ref).
"""
struct BSpline{D<:Degree} <: InterpolationType
    degree::D
end

bsplinetype(::Type{BSpline{D}}) where {D<:Degree} = D
bsplinetype(::BS) where {BS<:BSpline} = bsplinetype(BS)

degree(mode::BSpline) = mode.degree
degree(::NoInterp) = NoInterp()

iscomplete(mode::BSpline) = iscomplete(degree(mode))
iscomplete(deg::DegreeBC) = _iscomplete(deg.bc.gt)
iscomplete(deg::Degree) = true
_iscomplete(::Nothing) = false
_iscomplete(::GridType) = true

function Base.show(io::IO, bs::BSpline)
    print(io, "BSpline(")
    show(io, degree(bs))
    print(io, ')')
end

function Base.show(io::IO, deg::DegreeBC)
    print(io, nameof(typeof(deg)), '(')
    show(io, deg.bc)
    print(io, ')')
end

"""
    BSplineInterpolation{T,N,TCoefs,IT,Axs}

An interpolant-type for b-spline interpolation on a uniform grid with integer nodes.
`T` indicates the element type for operations like `collect(itp)`, and may also agree with
the values obtained from `itp(x, y, ...)` at least for certain types of `x` and `y`.
`N` is the dimensionality of the interpolant. The remaining type-parameters
describe the types of fields:
- the `coefs` field holds the interpolation coefficients. Depending on prefiltering,
  these may or may not be the same as the supplied array of interpolant values.
- `parentaxes` holds the axes of the parent. Depending on prefiltering this may be
  "narrower" than the axes of `coefs`.
- `it` holds the interpolation type, e.g., `BSpline(Linear())` or `(BSpline(Quadratic(OnCell()),BSpline(Linear()))`.

`BSplineInterpolation` objects are typically created with [`interpolate`](@ref).
However, for customized control you may also construct them with

    BSplineInterpolation(TWeights, coefs, it, axs)

where `T` gets computed from the product of `TWeights` and `eltype(coefs)`.
(This is equivalent to indicating that you'll be evaluating at locations `itp(x::TWeights, y::TWeights, ...)`.)
"""
struct BSplineInterpolation{T,N,TCoefs<:AbstractArray,IT<:DimSpec{BSpline},Axs<:Tuple{Vararg{AbstractUnitRange,N}}} <: AbstractInterpolation{T,N,IT}
    coefs::TCoefs
    parentaxes::Axs
    it::IT
end

function Base.:(==)(o1::BSplineInterpolation, o2::BSplineInterpolation)
    o1.it == o2.it &&
    o1.parentaxes == o2.parentaxes &&
    o1.coefs == o2.coefs
end

function BSplineInterpolation(::Type{TWeights}, A::AbstractArray{Tel,N}, it::IT, axs) where {N,Tel,TWeights<:Real,IT<:DimSpec{BSpline}}
    # String interpolation causes allocation, noinline avoids that unless they get called
    @noinline err_concrete(IT) = error("The b-spline type must be a concrete type (was $IT)")
    @noinline warn_concrete(A) = @warn("For performance reasons, consider using an array of a concrete type (typeof(A) == $(typeof(A)))")
    @noinline err_incomplete(it) = error("OnGrid/OnCell is not supplied for some of the interpolation modes in $it")
    @noinline err_singleton(A, it) = throw(ArgumentError("size $(size(A)) is inconsistent with $it, use NoInterp along singleton dimensions"))

    isconcretetype(IT) || err_concrete(IT)
    isconcretetype(typeof(A)) || warn_concrete(A)
    iscomplete(it) || err_incomplete(it)

    # Check that dimensions of size 1 are NoInterp
    is_singleton_ok(A, it) || err_singleton(A, it)

    # Compute the output element type when positions have type TWeights
    if isempty(A)
        T = Base.promote_op(*, TWeights, eltype(A))
    else
        T = typeof(zero(TWeights) * first(A))
    end
    BSplineInterpolation{T,N,typeof(A),IT,typeof(axs)}(A, fix_axis.(axs), it)
end

function BSplineInterpolation(A::AbstractArray{Tel,N}, it::IT, axs) where {N,Tel,IT<:DimSpec{BSpline}}
    @noinline err_axes(A, it, axs) = throw(ArgumentError("parentaxes $axs are inconsistent with coefficient array axes $(axes(A)) and interpolation type $it"))

    paxs = padded_axes(axs, it)
    all(map(⊆, paxs, axes(A))) || err_axes(A, it, axs)
    BSplineInterpolation(tweight(A), A, it, axs)
end

iscomplete(its::Tuple) = all(iscomplete, its)

is_singleton_ok(A, it) = is_singleton_ok(size(A), it)
@inline function is_singleton_ok(sz::Tuple, it)
    if sz[1] == 1
        it1 = getfirst(it)
        it1 isa NoInterp || return false
    end
    return is_singleton_ok(Base.tail(sz), getrest(it))
end
is_singleton_ok(::Tuple{}, it) = true

coefficients(itp::BSplineInterpolation) = itp.coefs
interpdegree(itp::BSplineInterpolation) = interpdegree(itpflag(itp))
interpdegree(::BSpline{T}) where T = T()
interpdegree(it::Tuple{Vararg{Union{BSpline,NoInterp},N}}) where N = interpdegree.(it)
itpflag(itp::BSplineInterpolation) = itp.it

size(itp::BSplineInterpolation) = map(length, itp.parentaxes)
axes(itp::BSplineInterpolation) = itp.parentaxes

lbounds(itp::BSplineInterpolation) = _lbounds(itp.parentaxes, itpflag(itp))
ubounds(itp::BSplineInterpolation) = _ubounds(itp.parentaxes, itpflag(itp))
_lbounds(axs, itp) = (lbound(axs[1], getfirst(itp)), _lbounds(Base.tail(axs), getrest(itp))...)
_ubounds(axs, itp) = (ubound(axs[1], getfirst(itp)), _ubounds(Base.tail(axs), getrest(itp))...)
_lbounds(::Tuple{}, itp) = ()
_ubounds(::Tuple{}, itp) = ()

lbound(ax::AbstractRange, bs::BSpline)   = lbound(ax, degree(bs))
lbound(ax::AbstractRange, deg::Degree)   = first(ax)
lbound(ax::AbstractRange, deg::DegreeBC) = lbound(ax, deg, deg.bc.gt)
ubound(ax::AbstractRange, bs::BSpline)   = ubound(ax, degree(bs))
ubound(ax::AbstractRange, deg::Degree)   = last(ax)
ubound(ax::AbstractRange, deg::DegreeBC) = ubound(ax, deg, deg.bc.gt)

lbound(ax::AbstractUnitRange, ::DegreeBC, ::OnCell) = first(ax) - 0.5
ubound(ax::AbstractUnitRange, ::DegreeBC, ::OnCell) = last(ax) + 0.5
lbound(ax::AbstractUnitRange, ::DegreeBC, ::OnGrid) = first(ax)
ubound(ax::AbstractUnitRange, ::DegreeBC, ::OnGrid) = last(ax)

fix_axis(r::Base.OneTo) = r
fix_axis(r::Base.Slice) = r
fix_axis(r::UnitRange) = Base.Slice(r)
fix_axis(r::AbstractUnitRange) = fix_axis(UnitRange(r))

count_interp_dims(::Type{BSI}, n) where BSI<:BSplineInterpolation = count_interp_dims(itptype(BSI), n)

function interpolate(::Type{TWeights}, ::Type{TC}, A, it::IT) where {TWeights,TC,IT<:DimSpec{BSpline}}
    Apad = prefilter(TWeights, TC, A, it)
    BSplineInterpolation(TWeights, Apad, it, axes(A))
end

"""
    itp = interpolate(A, interpmode)

Interpolate an array `A` in the mode determined by `interpmode`.
`interpmode` may be one of

- `NoInterp()`
- `BSpline(Constant())`
- `BSpline(Linear())`
- `BSpline(Quadratic(bc))` (see [`BoundaryCondition`](@ref))
- `BSpline(Cubic(bc))`

It may also be a tuple of such values, if you want to use different interpolation schemes along each axis.
"""
function interpolate(A::AbstractArray, it::IT) where {IT<:DimSpec{BSpline}}
    interpolate(tweight(A), tcoef(A), A, it)
end

# We can't just return a tuple-of-types due to julia #12500
tweight(A::AbstractArray) = Float64
tweight(A::AbstractArray{T}) where T<:AbstractFloat = T
tweight(A::AbstractArray{Rational{Int}}) = Rational{Int}
tweight(A::AbstractArray{T}) where {T<:Integer} = typeof(float(zero(T)))

tcoef(A::AbstractArray) = eltype(A)
tcoef(A::AbstractArray{Float32}) = Float32
tcoef(A::AbstractArray{Rational{Int}}) = Rational{Int}
tcoef(A::AbstractArray{T}) where {T<:Integer} = typeof(float(zero(T)))

function interpolate!(::Type{TWeights}, A::AbstractArray, it::IT) where {TWeights,IT<:DimSpec{BSpline}}
    # Set the bounds of the interpolant inward, if necessary
    axsA = axes(A)
    axspad = padded_axes(axsA, it)
    BSplineInterpolation(TWeights, prefilter!(TWeights, A, it), it, fix_axis.(padinset.(axsA, axspad)))
end
function interpolate!(A::AbstractArray, it::IT) where {IT<:DimSpec{BSpline}}
    interpolate!(tweight(A), A, it)
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
