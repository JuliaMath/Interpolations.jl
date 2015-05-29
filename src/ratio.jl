# Base's Rational type is too slow because of its penchant for calling gcd and div.
# So we roll our own.

immutable Ratio{T<:Integer} <: Real
    num::T
    den::T
end

convert{S}(::Type{BigFloat}, r::Ratio{S}) = BigFloat(r.num)/r.den
function convert{T<:FloatingPoint,S}(::Type{T}, r::Ratio{S})
    P = promote_type(T,S)
    convert(T, convert(P, r.num)/convert(P, r.den))
end
convert{T<:Integer}(::Type{Ratio{T}}, i::Integer) = Ratio{T}(convert(T, i), one(T))

*(x::Ratio, y::Ratio) = Ratio(x.num*y.num, x.den*y.den)
*(x::Ratio, y::Bool) = Ratio(x.num*y, x.den)
*(x::Ratio, y::Integer) = Ratio(x.num*y, x.den)
/(x::Ratio, y::Ratio) = Ratio(x.num*y.den, x.den*y.num)
/(x::Ratio, y::Integer) = Ratio(x.num, x.den*y)
+(x::Integer, y::Ratio) = Ratio(x*y.den + y.num, y.den)
-(x::Integer, y::Ratio) = Ratio(x*y.den - y.num, y.den)
+(x::Ratio, y::Ratio) = Ratio(x.num*y.den + x.den*y.num, x.den*y.den)
-(x::Ratio, y::Ratio) = Ratio(x.num*y.den - x.den*y.num, x.den*y.den)

promote_rule{T<:Integer,S<:Integer}(::Type{Ratio{T}}, ::Type{S}) = Ratio{promote_type(T,S)}
promote_rule{T<:Integer,S<:Integer}(::Type{Ratio{T}}, ::Type{Ratio{S}}) = Ratio{promote_type(T,S)}
promote_rule{T<:Integer,S<:FloatingPoint}(::Type{Ratio{T}}, ::Type{S}) = promote_type(T,S)

sqr(x) = x*x

