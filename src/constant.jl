type ConstantDegree <: Degree{0} end
type Constant{GR<:GridRepresentation} <: InterpolationType{ConstantDegree,None,GR} end
Constant{GR<:GridRepresentation}(::GR) = Constant{GR}()

function define_indices(::Constant, N)
    :(@nexprs $N d->(ix_d = clamp(round(Int,x_d), 1, size(itp,d))))
end

function coefficients(::Constant, N)
    :(@nexprs $N d->(c_d = one(typeof(x_d))))
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
