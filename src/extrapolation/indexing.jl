@generated function getindex{T}(exp::AbstractExtrapolation{T,1}, x)
    quote
        $(extrap_prep(exp, x))
        exp.itp[x]
    end
end

@generated function getindex{T,N,ITP,GT}(exp::AbstractExtrapolation{T,N,ITP,GT}, xs...)
    quote
        $(extrap_prep(exp, xs...))
        exp.itp[xs...]
    end
end
