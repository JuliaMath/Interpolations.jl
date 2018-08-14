module NonNumeric

# Test interpolation with a multi-valued type

using Interpolations, Test

import Base: +, -, *, /, ≈

struct MyPair{T}
    first::T
    second::T
end

# Here's the interface your type must define
(+)(p1::MyPair, p2::MyPair) = MyPair(p1.first+p2.first, p1.second+p2.second)
(-)(p1::MyPair, p2::MyPair) = MyPair(p1.first-p2.first, p1.second-p2.second)
(*)(n::Number, p::MyPair) = MyPair(n*p.first, n*p.second)
(*)(p::MyPair, n::Number) = n*p
(/)(p::MyPair, n::Number) = MyPair(p.first/n, p.second/n)
Base.zero(::Type{MyPair{T}}) where {T} = MyPair(zero(T),zero(T))
Base.promote_rule(::Type{MyPair{T1}}, ::Type{T2}) where {T1,T2<:Number} = MyPair{promote_type(T1,T2)}
≈(p1::MyPair, p2::MyPair) = (p1.first ≈ p2.first) & (p1.second ≈ p2.second)
# Base.promote_op(::typeof(*), ::Type{MyPair{T1}}, ::Type{T2}) where {T1,T2<:Number} = MyPair{promote_type(T1,T2)}
# Base.promote_op(::typeof(*), ::Type{T1}, ::Type{MyPair{T2}}) where {T1<:Number,T2} = MyPair{promote_type(T1,T2)}

# 1d
A0 = rand(20)
A = reinterpret(MyPair{Float64}, A0)
a1, a2 = A0[1:2:end], A0[2:2:end]
@test length(A) == 10
itp = interpolate(A, BSpline(Constant()), OnGrid())
@test itp(3.2) ≈ MyPair(A0[5],A0[6])
itp = interpolate(A, BSpline(Linear()), OnGrid())
@test itp(3.2) ≈ 0.8*MyPair(A0[5],A0[6]) + 0.2*MyPair(A0[7],A0[8])
it, gt = BSpline(Quadratic(Flat())), OnGrid()
itp = interpolate(A, it, gt)
@test itp(3.2) ≈ MyPair(interpolate(a1, it, gt)(3.2), interpolate(a2, it, gt)(3.2))

# 2d
A0 = rand(100)
A = reshape(reinterpret(MyPair{Float64}, A0), (10,5))
a1, a2 = reshape(A0[1:2:end], (10,5)), reshape(A0[2:2:end], (10,5))
for (it, gt) in ((BSpline(Constant()), OnGrid()),
                 (BSpline(Linear()), OnGrid()),
                 (BSpline(Quadratic(Flat())), OnGrid()))
    itp = interpolate(A, it, gt)
    @test itp(3.2,1.8) ≈ MyPair(interpolate(a1, it, gt)(3.2,1.8), interpolate(a2, it, gt)(3.2,1.8))
end

end
