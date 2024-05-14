export Lanczos

abstract type AbstractLanczos <: InterpolationType end

"""
    Lanczos{N}(a=4)

Lanczos resampling via a kernel with scale parameter `a` and support over `N` neighbors.

This form of interpolation is merely the discrete convolution of the samples
with a Lanczos kernel of size `a`. The size is directly related to how "far" the
interpolation will reach for information, and has `O(N^2)` impact on runtime.
An alternative implementation matching `lanczos4` from OpenCV is available as
Lanczos4OpenCV.
"""
struct Lanczos{N} <: AbstractLanczos
    a::Int

    function Lanczos{N}(a) where N
        N < a && @warn "Using a smaller support than scale for Lanczos window. Proceed with caution."
        new{N}(a)
    end
end

Lanczos(a=4) = Lanczos{a}(a)

"""
    LanczosInterpolation
"""
struct LanczosInterpolation{T,N,IT <: DimSpec{AbstractLanczos},A <: AbstractArray{T,N},P <: Tuple{Vararg{AbstractArray,N}}} <: AbstractInterpolation{T,N,IT}
    coefs::A
    parentaxes::P
    it::IT
end

LanczosInterpolation{T,N}(coefs::AbstractArray{T,N}, parentaxes::NTuple{N,AbstractArray}, it::IT) where {T,N,IT} =
    LanczosInterpolation{T,N,IT,typeof(coefs),typeof(parentaxes)}(coefs, parentaxes, it)

@generated degree(::Lanczos{N}) where {N} = :($N)

getknots(itp::LanczosInterpolation) = axes(itp)
coefficients(itp::LanczosInterpolation) = itp.coefs
itpflag(itp::LanczosInterpolation) = itp.it

size(itp::LanczosInterpolation) = map(length, itp.parentaxes)
axes(itp::LanczosInterpolation) = itp.parentaxes
lbounds(itp::LanczosInterpolation) = map(first, itp.parentaxes)
ubounds(itp::LanczosInterpolation) = map(last, itp.parentaxes)

function interpolate(A::AbstractArray{T}, it::AbstractLanczos) where T
    Apad = copy_with_padding(float(T), A, it)
    return LanczosInterpolation(Apad, axes(A), it)
end

@inline function (itp::LanczosInterpolation{T,N})(x::Vararg{Number,N}) where {T,N}
    @boundscheck (checkbounds(Bool, itp, x...) || Base.throw_boundserror(itp, x))
    wis = weightedindexes((value_weights,), itpinfo(itp)..., x)
    InterpGetindex(itp)[wis...]
end

function weightedindex_parts(fs, it::AbstractLanczos, ax::AbstractUnitRange{<:Integer}, x)
    pos, δx = positions(it, ax, x)
    (position = pos, coefs = fmap(fs, it, δx))
end

function positions(it::AbstractLanczos, ax, x)
    xf = floorbounds(x, ax)
    δx = x - xf
    fast_trunc(Int, xf) - degree(it) + 1, δx
end

function value_weights(it::Lanczos, δx::S) where S
    N = degree(it)
    # short-circuit if integral
    isinteger(δx) && return ntuple(i->convert(float(S), i == N - δx), Val(2N))

    # LUTs
    #it.a === N === 4 && return _lanczos4(δx)

    cs = ntuple(i -> lanczos(N - i + δx, it.a, N), Val(2N))
    sum_cs = sum(cs)
    normed_cs = ntuple(i -> cs[i] / sum_cs, Val(length(cs)))
    return normed_cs
end

function padded_axis(ax::AbstractUnitRange, it::AbstractLanczos)
    N = degree(it)
    return first(ax) - N + 1:last(ax) + N
end

# precise implementations for fast evaluation of common kernels

"""
    lanczos(x, a, n=a)

Implementation of the [Lanczos kernel](https://en.wikipedia.org/wiki/Lanczos_resampling)
"""
lanczos(x::T, a::Integer, n=a) where {T} = abs(x) < n ? T(sinc(x) * sinc(x / a)) : zero(T)
