using Base.Cartesian

export Lanczos

"""
    Lanczos(a=4)

Lanczos resampling via a kernel with size `a`.

This form of interpolation is merely the discrete convolution of the samples with a Lanczos kernel of size `a`. The size is directly related to how "far" the interpolation will reach for information. A default value of 4 matches the OpenCV implementation `lanczos4`.
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

struct LanczosInterpolation{T, N, IT<:DimSpec{Lanczos}} <: AbstractInterpolation{T, N, IT}
    data::Array{T,N}
    degree::Int
    it::IT
end

@generated function (lz::LanczosInterpolation{T, N})(pts...) where {T,N}
    quote
        a = lz.degree
        K = lz.data
        w = s = 0
        # get weighting
        @nloops $N i (d -> -a+1:a+1) begin
            q = 1
            ls = @nexprs $N d -> q *= lanczos(T, i_d - pts[d] + floor(pts[d]), a)
            w += q
        end
        # get value
        @nloops $N i (d -> -a+1:a+1) (d -> idx_d = floor(Int, pts[d]) + i_d) begin
            (@ncall $N checkbounds Bool K d -> idx_d) || continue
            f = @nref $N K idx
            ls = @nexprs $N d -> f *= lanczos(T, i_d - pts[d] + floor(pts[d]), a)

            s += f
        end
        return s / w
    end
end

getknots(lz::LanczosInterpolation) = axes(lz.data)

size(lz::LanczosInterpolation) = size(lz.data)
size(lz::LanczosInterpolation, i) = size(lz.data, i)
axes(lz::LanczosInterpolation) = axes(lz.data)
axes(lz::LanczosInterpolation, i) = axes(lz.data, i)
lbounds(lz::LanczosInterpolation) = map(first, axes(lz.data))
lbounds(lz::LanczosInterpolation, i) = first(axes(lz.data, i))
ubounds(lz::LanczosInterpolation) = map(last, axes(lz.data))
ubounds(lz::LanczosInterpolation, i) = last(axes(lz.data, i))
