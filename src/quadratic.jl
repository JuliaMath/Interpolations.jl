type QuadraticDegree <: Degree{2} end
type Quadratic{BC<:BoundaryCondition,GR<:GridRepresentation} <: InterpolationType{QuadraticDegree,BC,GR} end

Quadratic{BC<:BoundaryCondition,GR<:GridRepresentation}(::BC, ::GR) = Quadratic{BC,GR}()

function bc_gen(q::Quadratic, N)
    quote
        pad = padding($q)
        @nexprs $N d->(ix_d = clamp(@compat round(Integer, x_d), 1, size(itp,d)) + pad)
    end
end

function indices(q::Quadratic, N)
    quote
        pad = padding($q)
        @nexprs $N d->begin
            ixp_d = ix_d + 1
            ixm_d = ix_d - 1

            fx_d = x_d - convert(typeof(x_d), ix_d - pad)
        end
    end
end

function coefficients(::Quadratic, N)
    quote
        @nexprs $N d->begin
            cm_d = .5 * (fx_d-.5)^2
            c_d = .75 - fx_d^2
            cp_d = .5 * (fx_d+.5)^2
        end
    end
end

# This assumes integral values ixm_d, ix_d, and ixp_d,
# coefficients cm_d, c_d, and cp_d, and an array itp.coefs
function index_gen(degree::QuadraticDegree, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        symm, sym, symp =  symbol(string("cm_",d)), symbol(string("c_",d)), symbol(string("cp_",d))
        return :($symm * $(index_gen(degree, N, offsets...,-1)) + $sym * $(index_gen(degree, N, offsets..., 0)) +
                 $symp * $(index_gen(degree, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

# Quadratic interpolation has 1 extra coefficient at each boundary
# Therefore, make the coefficient array 1 larger in each direction,
# in each dimension.
padding(::Quadratic) = 1

function inner_system_diags{T}(::Type{T}, n::Int, ::Quadratic)
    du = fill(convert(T,1/8), n-1)
    d = fill(convert(T,3/4),n)
    dl = copy(du)
    d[1] = d[end] = du[1] = dl[end] = convert(T,0)
    (dl,d,du)
end

function prefiltering_system{T}(::Type{T}, n::Int, q::Quadratic{Flat,OnCell})
    dl,d,du = inner_system_diags(T,n,q)
    d[1] = d[end] = -1
    du[1] = dl[end] = 1
    lufact!(Tridiagonal(dl, d, du)), zeros(T, n)
end

function prefiltering_system{T}(::Type{T}, n::Int, q::Quadratic{Flat,OnGrid})
    dl,d,du = inner_system_diags(T,n,q)
    d[1] = d[end] = convert(T, -1)
    du[1] = dl[end] = convert(T, 0)

    rowspec = zeros(T,n,2)
    # first row     last row
    rowspec[1,1] = rowspec[n,2] = convert(T, 1)
    colspec = zeros(T,2,n)
    # third col     third-to-last col
    colspec[1,3] = colspec[2,n-2] = convert(T, 1)
    valspec = zeros(T,2,2)
    # val for [1,3], val for [n,n-2]
    valspec[1,1] = valspec[2,2] = convert(T, 1)

    Woodbury(lufact!(Tridiagonal(dl, d, du)), rowspec, valspec, colspec), zeros(T, n)
end

function prefiltering_system{T}(::Type{T}, n::Int, q::Quadratic{LinearBC})
    dl,d,du = inner_system_diags(T,n,q)
    d[1] = d[end] = convert(T, 1)
    du[1] = dl[end] = convert(T, -2)

    rowspec = zeros(T,n,2)
    # first row     last row
    rowspec[1,1] = rowspec[n,2] = convert(T, 1)
    colspec = zeros(T,2,n)
    # third col     third-to-last col
    colspec[1,3] = colspec[2,n-2] = convert(T, 1)
    valspec = zeros(T,2,2)
    # val for [1,3], val for [n,n-2]
    valspec[1,1] = valspec[2,2] = convert(T, 1)

    Woodbury(lufact!(Tridiagonal(dl, d, du)), rowspec, valspec, colspec), zeros(T, n)
end
