using Base.Cartesian
using StaticArrays

export Lanczos

"""
    Lanczos(a=4, n=a)

Lanczos resampling via a kernel with scale parameter `a` and support over `n` neighbors.

This form of interpolation is merely the discrete convolution of the samples with a Lanczos kernel of size `a`. The size is directly related to how "far" the interpolation will reach for information, and has `O(n^2)` impact on runtime. A default value of 4 matches the OpenCV implementation `lanczos4`.
"""
struct Lanczos <: InterpolationType
    a::Int
    n::Int
    
    function Lanczos(a, n)
        n < a && @warn "Using a smaller support than scale for Lanczos window. Proceed with caution."
        new(a, n)
    end
end

Lanczos(a=4) = Lanczos(a, a)


"""
    lanczos(x, a, n=a)

Implementation of the [Lanczos kernel](https://en.wikipedia.org/wiki/Lanczos_resampling)
"""
lanczos(x::T, a::Integer, n=a) where {T} = abs(x) < n ? T(sinc(x) * sinc(x / a)) : zero(T)

struct LanczosInterpolation{T,N,IT <: DimSpec{Lanczos},A <: AbstractArray{T,N},P <: Tuple{Vararg{AbstractArray,N}}} <: AbstractInterpolation{T,N,IT}
    coefs::A
    parentaxes::P
    it::IT
end

getknots(itp::LanczosInterpolation) = axes(itp)
coefficients(itp::LanczosInterpolation) = itp.coefs
itpflag(itp::LanczosInterpolation) = itp.it

size(itp::LanczosInterpolation) = map(length, itp.parentaxes)
axes(itp::LanczosInterpolation) = itp.parentaxes
lbounds(itp::LanczosInterpolation) = map(first, itp.parentaxes)
ubounds(itp::LanczosInterpolation) = map(last, itp.parentaxes)

function interpolate(A::AbstractArray{T}, it::Lanczos) where T
    Apad = copy_with_padding(T, A, it)
    return LanczosInterpolation(Apad, axes(A), it)
end

@inline function (itp::LanczosInterpolation{T,N})(x::Vararg{<:Number,N}) where {T,N}
    @boundscheck (checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x))
    wis = weightedindexes((value_weights,), itpinfo(itp)..., x)
    coefficients(itp)[wis...]
end

function weightedindex_parts(fs::F, it::Lanczos, ax::AbstractUnitRange{<:Integer}, x) where F
    pos, δx = positions(it, ax,  x)
    (position = pos, coefs = fmap(fs, it, δx))
end

function positions(it::Lanczos, ax, x)
    xf = floorbounds(x, ax)
    δx = x - xf
    fast_trunc(Int, xf) - it.n + 1, δx
end

function value_weights(it::Lanczos, δx::T) where T
    # LUTs
    it.a == it.n == 4 && return _lanczos4(δx)

    # short-circuit if integral
    isinteger(δx) && return ntuple(i -> i == it.n - δx ? one(T) : zero(T), 2it.n)

    idxs = -it.n + 1:it.n
    cs = @. lanczos(idxs + δx, it.a, it.n)
    return Tuple(cs ./ sum(cs))
end

function padded_axis(ax::AbstractUnitRange, it::Lanczos)
    return first(ax) - it.n + 1:last(ax) + it.n
end

# precise implementations for fast evaluation of common kernels

const s45 = 0.70710678118654752440084436210485
const l4_2d_cs = SA[1 0; -s45 -s45; 0 1; s45 -s45; -1 0; s45 s45; 0 -1; -s45 s45]

function _lanczos4(δx::T) where T
    isinteger(δx) && return ntuple(i -> i == 4 - δx ? one(T) : zero(T), 8)
    y0 = -(δx + 3) * π/4
    s0, c0 = sincos(y0)
    cs = ntuple(8) do i
        y = -(δx + 3 - i) * π / 4
        (l4_2d_cs[i, 1] * s0 + l4_2d_cs[i, 2] * c0) / y^2
    end
    return cs ./ sum(cs)
end
