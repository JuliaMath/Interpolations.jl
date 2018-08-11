module NonNumeric

# Test interpolation with a multi-valued type

using Interpolations

import Base: +, -, *, /

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
Base.promote_op(::typeof(*), ::Type{MyPair{T1}}, ::Type{T2}) where {T1,T2<:Number} = MyPair{promote_type(T1,T2)}
Base.promote_op(::typeof(*), ::Type{T1}, ::Type{MyPair{T2}}) where {T1<:Number,T2} = MyPair{promote_type(T1,T2)}

# 1d
A = reinterpret(MyPair{Float64}, rand(20))
itp = interpolate(A, BSpline(Constant()), OnGrid())
itp[3.2]
itp = interpolate(A, BSpline(Linear()), OnGrid())
itp[3.2]
itp = interpolate(A, BSpline(Quadratic(Flat())), OnGrid())
itp[3.2]

# 2d
A = reshape(reinterpret(MyPair{Float64}, rand(100)), (10,5))
itp = interpolate(A, BSpline(Constant()), OnGrid())
itp[3.2,1.8]
itp = interpolate(A, BSpline(Linear()), OnGrid())
itp[3.2,1.8]
itp = interpolate(A, BSpline(Quadratic(Flat())), OnGrid())
itp[3.2,1.8]

end
