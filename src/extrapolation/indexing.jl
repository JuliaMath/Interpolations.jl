@generated function getindex{T}(etp::AbstractExtrapolation{T,1}, x)
    quote
        $(extrap_prep(etp, x))
        etp.itp[x]
    end
end

@generated function getindex{T,N,ITP,GT}(etp::AbstractExtrapolation{T,N,ITP,GT}, xs...)
    quote
        $(extrap_prep(etp, xs...))
        etp.itp[xs...]
    end
end
