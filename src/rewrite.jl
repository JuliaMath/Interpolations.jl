struct MyInterp{T,N,A<:AbstractArray{T,N}}
    data::A
end

@inline (itp::MyInterp{T,N})(i::Vararg{<:Number,N}) where {T,N} = expand(itp, i)

@inline function expand(itp::MyInterp, ipre::Tuple{Vararg{Number,L}}, ipost::Vararg{Integer,M}) where {L,M}  # force specialization
    ifront, ilast = Base.front(ipre), ipre[end]
    im, ip = floor(ilast), ceil(ilast)
    return (ip - ilast)*expand(itp, ifront, unsafe_trunc(Int, im), ipost...) +
           (ilast - im)*expand(itp, ifront, unsafe_trunc(Int, ip), ipost...)
end

@inline expand(itp::MyInterp, ::Tuple{}, ipost::Vararg{Integer,N}) where N =
    @inbounds itp.data[CartesianIndex(ipost)]
