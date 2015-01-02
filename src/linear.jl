type LinearDegree <: Degree{1} end
type Linear{GR<:GridRepresentation} <: InterpolationType{LinearDegree,None,GR} end
Linear{GR<:GridRepresentation}(::GR) = Linear{GR}()

function define_indices(::Linear, N)
    quote
        @nexprs $N d->begin
            ix_d = clamp(floor(Int,x_d), 1, size(itp,d)-1)
            ixp_d = ix_d + 1
            fx_d = x_d - convert(typeof(x_d), ix_d)
        end
    end
end

function coefficients(l::Linear, N)
    :(@nexprs $N d->($(coefficients(l, N, :d))))
end

function coefficients(::Linear, N, d)
    sym, symp, symfx = symbol(string("c_",d)), symbol(string("cp_",d)), symbol(string("fx_",d))
    quote
        $sym = one(typeof($symfx)) - $symfx
        $symp = $symfx
    end
end

function gradient_coefficients(::Linear,N,d)
    sym, symp, symfx = symbol(string("c_",d)), symbol(string("cp_",d)), symbol(string("fx_",d))
    quote
        $sym = -one(typeof($symfx))
        $symp = one(typeof($symfx))
    end
end

# This assumes fractional values 0 <= fx_d <= 1, integral values ix_d and ixp_d (typically ixp_d = ix_d+1,
#except at boundaries), and an array itp.coefs
function index_gen(degree::LinearDegree, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        sym = symbol("c_"*string(d))
        symp = symbol("cp_"*string(d))
        return :($sym * $(index_gen(degree, N, offsets..., 0)) + $symp * $(index_gen(degree, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
