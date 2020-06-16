using Base.Cartesian

export Lanczos

"""
    Lanczos(a=4)

Lanczos resampling via a kernel with size `a`.

This form of interpolation is merely the discrete convolution of the samples with a Lanczos kernel of size `a`. The size is directly related to how "far" the interpolation will reach for information, and has `O(a^2)` impact on runtime. A default value of 4 matches the OpenCV implementation `lanczos4`.
"""
struct Lanczos <: InterpolationType
    a::Integer
end

Lanczos() = Lanczos(4)

"""
    lanczos(x, a)

Implementation of the [Lanczos kernel](https://en.wikipedia.org/wiki/Lanczos_resampling)
"""
lanczos(x::T, a::Integer) where{T} = lanczos(T, x, a)
lanczos(T, x, a::Integer) = abs(x) < a ? T(sinc(x) * sinc(x / a)) : zero(T)

interpolate(s::AbstractArray, it::Lanczos) = LanczosInterpolation(s, it.a, it)

struct LanczosInterpolation{T, N, IT<:DimSpec{Lanczos},A<:AbstractArray{T, N}} <: AbstractInterpolation{T, N, IT}
    data::A
    degree::Int
    it::IT
end

@inline @generated function (itp::LanczosInterpolation{T, N})(pts::Vararg{<:Number, N}) where {T,N}
    quote
        # if we can just index, let's
        all(isinteger, pts) && return K[convert.(Int, pts)...]
        @boundscheck checkbounds(K, pts...)
        
        a = itp.degree
        K = itp.data
        w = s = 0

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

getknots(itp::LanczosInterpolation) = axes(itp.data)
coefficients(itp::LanczosInterpolation) = itp.data

size(itp::LanczosInterpolation) = size(itp.data)
size(itp::LanczosInterpolation, i) = size(itp.data, i)
axes(itp::LanczosInterpolation) = axes(itp.data)
axes(itp::LanczosInterpolation, i) = axes(itp.data, i)
lbounds(itp::LanczosInterpolation) = map(first, axes(itp.data))
lbounds(itp::LanczosInterpolation, i) = first(axes(itp.data, i))
ubounds(itp::LanczosInterpolation) = map(last, axes(itp.data))
ubounds(itp::LanczosInterpolation, i) = last(axes(itp.data, i))
