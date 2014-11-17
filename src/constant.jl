type ConstantDegree <: Degree{0} end
type Constant{GR<:GridRepresentation} <: InterpolationType{ConstantDegree,BC.None,GR} end

Constant{GR<:GridRepresentation}(::GR) = Constant{GR}()

function bc_gen{IT<:Constant}(::IT, N)
    quote
        @nexprs $N d->(ix_d = iround(x_d))
    end
end

function indices(::ConstantDegree, N)
    quote
        # Constant interpolation doesn't need an fx_d
    end
end

function coefficients(::ConstantDegree, N)
    quote
        @nexprs $N d->(c_d = one(typeof(x_d)))
    end
end

function index_gen(degree::ConstantDegree, N::Integer, offsets...)
    if (length(offsets) < N)
        d = length(offsets)+1
        sym = symbol("c_"*string(d))
        return :($sym * $(index_gen(degree, N, offsets..., 0)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end
