using Base.Cartesian
using StaticArrays

export Lanczos

"""
    Lanczos(a)
    Lanczos{a}()

Lanczos resampling via a kernel with size `a`.

This form of interpolation is merely the discrete convolution of the samples with a Lanczos kernel of size `a`. The size is directly related to how "far" the interpolation will reach for information, and has `O(a^2)` impact on runtime. A default value of 4 matches the OpenCV implementation `lanczos4`.
"""
struct Lanczos{A} <: InterpolationType end

Lanczos(a=4) = Lanczos{a}()

"""
    lanczos(x, a)

Implementation of the [Lanczos kernel](https://en.wikipedia.org/wiki/Lanczos_resampling)
"""
lanczos(x::T, a::Integer) where{T} = lanczos(T, x, a)
lanczos(T, x, a::Integer) = abs(x) < a ? T(sinc(x) * sinc(x / a)) : zero(T)

@generated function value_weights(::Lanczos{A}, δx) where A
    quote
        idxs = -$A+1:$A
        cs = lanczos.(idxs .+ δx, $A)
        return tuple((cs ./ sum(cs))...)
    end
end
@generated degree(::Lanczos{A}) where {A} = :($A)

# @generated function


interpolate(s::AbstractArray, it::Lanczos) = LanczosInterpolation(axes(s), s, it)

struct LanczosInterpolation{T, N, IT<:DimSpec{Lanczos},A<:AbstractArray{T, N},P<:Tuple{Vararg{AbstractArray}}} <: AbstractInterpolation{T, N, IT}
    parentaxes::P
    coefs::A
    it::IT
end

# @generated degree(lz::LanczosInterpolation{T,N,IT{A}}) where {T,N,A,IT} = $A
getknots(itp::LanczosInterpolation) = itp.parentaxes
coefficients(itp::LanczosInterpolation) = itp.coefs

size(itp::LanczosInterpolation) = size(itp.coefs)
size(itp::LanczosInterpolation, i) = size(itp.coefs, i)
axes(itp::LanczosInterpolation) = axes(itp.coefs)
axes(itp::LanczosInterpolation, i) = axes(itp.coefs, i)
lbounds(itp::LanczosInterpolation) = map(first, axes(itp.coefs))
lbounds(itp::LanczosInterpolation, i) = first(axes(itp.coefs, i))
ubounds(itp::LanczosInterpolation) = map(last, axes(itp.coefs))
ubounds(itp::LanczosInterpolation, i) = last(axes(itp.coefs, i))



@inline @generated function (itp::LanczosInterpolation{T, N})(pts::Vararg{<:Number, N}) where {T,N}
    quote
        # if we can just index, let's
        all(isinteger, pts) && return itp.coefs[convert.(Int, pts)...]

        a = itp.degree
        K = itp.coefs
        w = s = 0
        @boundscheck checkbounds(K, pts...)

        @inbounds @nloops $N i (d -> -a+1:a+1) (d -> idx_d = floor(Int, pts[d]) + i_d) begin
            (@ncall $N checkbounds Bool K d -> idx_d) || continue
            # get weighting
            q = 1
            @nexprs $N d -> q *= lanczos(T, i_d - pts[d] + floor(pts[d]), a)
            w += q
            # get value
            s += (@nref $N K idx) * q
        end
        return s / w
    end
end
