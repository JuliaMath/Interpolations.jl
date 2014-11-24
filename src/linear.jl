type LinearDegree <: Degree{1} end
type Linear{GR<:GridRepresentation} <: InterpolationType{LinearDegree,None,GR} end
Linear{GR<:GridRepresentation}(::GR) = Linear{GR}()

Interpolation{T,N,EB<:ExtrapolationBehavior}(A::Array{T,N}, ::Linear, ::EB) = Interpolation{T,N,Linear{OnGrid},EB}(A)

function bc_gen(::Linear{OnGrid}, N)
    quote
        # LinearOnGrid's boundary condition doesn't require any action;
        # just return the cell index of the interpolation coordinate
        @nexprs $N d->(ix_d = ifloor(x_d))
    end
end

function indices(::Linear{OnGrid}, N)
    quote
        # fx_d is a parameter in [0,1] such that x_d = ix_d + fx_d
        @nexprs $N d->(fx_d = x_d - convert(typeof(x_d), ix_d))
        # ixp_d is the index in dimension d of the nearest node *after* the interpolation point
        @nexprs $N d->(ixp_d = ix_d + 1)
    end
end

function coefficients(::Linear{OnGrid}, N)
    quote
        @nexprs $N d->begin
            c_d = one(typeof(fx_d)) - fx_d
            cp_d = fx_d
        end
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
