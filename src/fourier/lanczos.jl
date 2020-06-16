using Base.Cartesian
using StaticArrays

export Lanczos

abstract type Kernel <: InterpolationType end

"""
    Lanczos(a=4)

Lanczos resampling via a kernel with size `a`.

This form of interpolation is merely the discrete convolution of the samples with a Lanczos kernel of size `a`. The size is directly related to how "far" the interpolation will reach for information, and has `O(a^2)` impact on runtime. A default value of 4 matches the OpenCV implementation `lanczos4`.
"""
struct Lanczos <: Kernel
    a::Integer

    # function Lanczos(a)
    #     a == 3 ? Lanczos3() :
    #     a == 4 ? Lanczos4() :
    #              new(a)
    # end
end

struct Lanczos4 <: Kernel end
struct Lanczos3 <: Kernel end

Lanczos() = Lanczos(4)

"""
    lanczos(x, a)

Implementation of the [Lanczos kernel](https://en.wikipedia.org/wiki/Lanczos_resampling)
"""
lanczos(x::T, a::Integer) where{T} = lanczos(T, x, a)
lanczos(T, x, a::Integer) = abs(x) < a ? T(sinc(x) * sinc(x / a)) : zero(T)

interpolate(s::AbstractArray, it::Lanczos) = LanczosInterpolation(s, it.a, it)

struct LanczosInterpolation{T, N, IT<:DimSpec{Kernel},A<:AbstractArray{T, N}} <: AbstractInterpolation{T, N, IT}
    data::A
    degree::Int
    it::IT
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

@inline @generated function (itp::LanczosInterpolation{T, N})(pts::Vararg{<:Number, N}) where {T,N}
    quote
        # if we can just index, let's
        all(isinteger, pts) && return itp.data[convert.(Int, pts)...]

        a = itp.degree
        K = itp.data
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

# const s45 = 0.70710678118654752440084436210485
# const l4_2d_coeffs = SArray[1    0  ;
#                            -s45 -s45;
#                             0    1  ;
#                             s45 -s45;
#                            -1    0  ;
#                             s45  s45;
#                             0   -1  ;
#                            -s45  s45]

# # specializations
# @inline function (itp::LanczosInterpolation{T,1,Lanczos4})(pts::Vararg{<:Number, 1}) where T
# end

# function _getcoeffs(x)
#     iszero(x) && return MArray[0, 0, 0, 1, 0, 0, 0, 0]
#     y0 = -(x + 3) * π / 4
#     s0, c0 = sincos(y0)
#     cs = MArray(undef, 8)
#     for i in eachindex(cs)
#         y = -(x + 3 - i) * π / 4
#         cs[i] = (l4_2d_coeffs[i, 1] * s0 + l4_2d_coeffs[i, 2] * c0) / y^2
#     end
#     s = 1 / sum(cs)
#     cs .*= s
#     return cs
# end

# const l4_2d_table = begin
#     scale = 1 / tabsz
    
# end
