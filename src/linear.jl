type LinearDegree <: Degree{1} end

type LinearOnGrid <: InterpolationType{LinearDegree,BCnone,OnGrid} end
typealias Linear LinearOnGrid

function interp_gen(::Type{Linear}, N)
    quote
        # ix_d is the index in dimension d of the nearest node *before* the interpolation point
        @nexprs $N d->(ix_d = ifloor(x_d))
        # fx_d is a parameter in [0,1] such that x_d = ix_d + fx_d
        @nexprs $N d->(fx_d = x_d - convert(typeof(x_d), ix_d))
        # ixp_d is the index in dimension d of the nearest node *after* the interpolation point
        @nexprs $N d->(ixp_d = ix_d + 1)

        # index_gen generates, recursively, an indexing expression that interpolates the data
        @inbounds ret = $(index_gen(Linear, N))
        ret
    end
end

#body_gen{BC<:BoundaryCondition,EB<:ExtrapolationBehavior}(::Type{Linear}, ::Type{BC}, ::Type{EB}, N) = body_gen(Linear, EB, N)

offsetsym(off, d) = off == -1 ? symbol(string("ixm_", d)) :
                    off ==  0 ? symbol(string("ix_", d)) :
                    off ==  1 ? symbol(string("ixp_", d)) :
                    off ==  2 ? symbol(string("ixpp_", d)) : error("offset $off not recognized")

# This assumes fractional values 0 <= fx_d <= 1, integral values ix_d and ixp_d (typically ixp_d = ix_d+1,
#except at boundaries), and an array itp.coefs
function index_gen(::Type{Linear}, N::Integer, offsets...)
    if length(offsets) < N
        sym = symbol("fx_"*string(length(offsets)+1))
        return :((one($sym)-$sym) * $(index_gen(Linear, N, offsets..., 0)) + $sym * $(index_gen(Linear, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
