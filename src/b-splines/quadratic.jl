immutable Quadratic{BC<:BoundaryCondition} <: Degree{2} end
Quadratic{BC<:BoundaryCondition}(::Type{BC}) = Quadratic{BC}

function define_indices_d{BC}(::Type{BSpline{Quadratic{BC}}}, d, pad)
    symix, symixm, symixp = symbol("ix_",d), symbol("ixm_",d), symbol("ixp_",d)
    symx, symfx = symbol("x_",d), symbol("fx_",d)
    quote
        # ensure that all three ix_d, ixm_d, and ixp_d are in-bounds no matter
        # the value of pad
        $symix = clamp(round(Int, real($symx)), 2-$pad, size(itp,$d)+$pad-1)
        $symfx = $symx - $symix
        $symix += $pad # padding for oob coefficient
        $symixp = $symix + 1
        $symixm = $symix - 1
    end
end
function define_indices_d(::Type{BSpline{Quadratic{Periodic}}}, d, pad)
    symix, symixm, symixp = symbol("ix_",d), symbol("ixm_",d), symbol("ixp_",d)
    symx, symfx = symbol("x_",d), symbol("fx_",d)
    quote
        $symix = clamp(round(Int, real($symx)), 1, size(itp,$d))
        $symfx = $symx - $symix
        $symixp = mod1($symix + 1, size(itp,$d))
        $symixm = mod1($symix - 1, size(itp,$d))
    end
end
function define_indices_d(::Type{BSpline{Quadratic{InPlace}}}, d, pad)
    symix, symixm, symixp = symbol("ix_",d), symbol("ixm_",d), symbol("ixp_",d)
    symx, symfx = symbol("x_",d), symbol("fx_",d)
    pad == 0 || error("Use InPlace only with interpolate!")
    quote
        # ensure that all three ix_d, ixm_d, and ixp_d are in-bounds no matter
        # the value of pad
        $symix = clamp(round(Int, real($symx)), 1, size(itp,$d))
        $symfx = $symx - $symix
        $symix += $pad # padding for oob coefficient
        $symixp = min(size(itp,$d), $symix + 1)
        $symixm = max(1, $symix - 1)
    end
end

function coefficients{Q<:Quadratic}(::Type{BSpline{Q}}, N, d)
    symm, sym, symp =  symbol(string("cm_",d)), symbol(string("c_",d)), symbol(string("cp_",d))
    symfx = symbol(string("fx_",d))
    quote
        $symm = sqr($symfx - SimpleRatio(1,2))/2
        $sym  = SimpleRatio(3,4) - sqr($symfx)
        $symp = sqr($symfx + SimpleRatio(1,2))/2
    end
end

function gradient_coefficients{Q<:Quadratic}(::Type{BSpline{Q}}, d)
    symm, sym, symp =  symbol(string("cm_",d)), symbol(string("c_",d)), symbol(string("cp_",d))
    symfx = symbol(string("fx_",d))
    quote
        $symm = $symfx - SimpleRatio(1,2)
        $sym = -2 * $symfx
        $symp = $symfx + SimpleRatio(1,2)
    end
end

# This assumes integral values ixm_d, ix_d, and ixp_d,
# coefficients cm_d, c_d, and cp_d, and an array itp.coefs
function index_gen{Q<:Quadratic,IT<:DimSpec{BSpline}}(::Type{BSpline{Q}}, ::Type{IT}, N::Integer, offsets...)
    if length(offsets) < N
        d = length(offsets)+1
        symm, sym, symp =  symbol(string("cm_",d)), symbol(string("c_",d)), symbol(string("cp_",d))
        return :($symm * $(index_gen(IT, N, offsets...,-1)) + $sym * $(index_gen(IT, N, offsets..., 0)) +
                 $symp * $(index_gen(IT, N, offsets..., 1)))
    else
        indices = [offsetsym(offsets[d], d) for d = 1:N]
        return :(itp.coefs[$(indices...)])
    end
end

padding{BC<:BoundaryCondition}(::Type{BSpline{Quadratic{BC}}}) = Val{1}()
padding(::Type{BSpline{Quadratic{Periodic}}}) = Val{0}()

function inner_system_diags{T,Q<:Quadratic}(::Type{T}, n::Int, ::Type{Q})
    du = fill(convert(T, SimpleRatio(1,8)), n-1)
    d = fill(convert(T, SimpleRatio(3,4)), n)
    dl = copy(du)
    (dl,d,du)
end

function prefiltering_system{T,TCoefs,BC<:Union(Flat,Reflect)}(::Type{T}, ::Type{TCoefs}, n::Int, ::Type{Quadratic{BC}}, ::Type{OnCell})
    dl,d,du = inner_system_diags(T,n,Quadratic{BC})
    d[1] = d[end] = -1
    du[1] = dl[end] = 1
    lufact!(Tridiagonal(dl, d, du), Val{false}), zeros(TCoefs, n)
end

function prefiltering_system{T,TCoefs}(::Type{T}, ::Type{TCoefs}, n::Int, ::Type{Quadratic{InPlace}}, ::Type{OnCell})
    dl,d,du = inner_system_diags(T,n,Quadratic{InPlace})
    d[1] = d[end] = convert(T, SimpleRatio(7,8))
    lufact!(Tridiagonal(dl, d, du), Val{false}), zeros(TCoefs, n)
end

function prefiltering_system{T,TCoefs,BC<:Union(Flat,Reflect)}(::Type{T}, ::Type{TCoefs}, n::Int, ::Type{Quadratic{BC}}, ::Type{OnGrid})
    dl,d,du = inner_system_diags(T,n,Quadratic{BC})
    d[1] = d[end] = -1
    du[1] = dl[end] = 0

    rowspec = spzeros(T,n,2)
    # first row     last row
    rowspec[1,1] = rowspec[n,2] = 1
    colspec = spzeros(T,2,n)
    # third col     third-to-last col
    colspec[1,3] = colspec[2,n-2] = 1
    valspec = zeros(T,2,2)
    # [1,3]         [n,n-2]
    valspec[1,1] = valspec[2,2] = 1

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), rowspec, valspec, colspec), zeros(TCoefs, n)
end

function prefiltering_system{T,TCoefs,GT<:GridType}(::Type{T}, ::Type{TCoefs}, n::Int, ::Type{Quadratic{Line}}, ::Type{GT})
    dl,d,du = inner_system_diags(T,n,Quadratic{Line})
    d[1] = d[end] = 1
    du[1] = dl[end] = -2

    rowspec = spzeros(T,n,2)
    # first row     last row
    rowspec[1,1] = rowspec[n,2] = 1
    colspec = spzeros(T,2,n)
    # third col     third-to-last col
    colspec[1,3] = colspec[2,n-2] = 1
    valspec = zeros(T,2,2)
    # [1,3]         [n,n-2]
    valspec[1,1] = valspec[2,2] = 1

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), rowspec, valspec, colspec), zeros(TCoefs, n)
end

function prefiltering_system{T,TCoefs,GT<:GridType}(::Type{T}, ::Type{TCoefs}, n::Int, ::Type{Quadratic{Free}}, ::Type{GT})
    dl,d,du = inner_system_diags(T,n,Quadratic{Free})
    d[1] = d[end] = 1
    du[1] = dl[end] = -3

    rowspec = spzeros(T,n,4)
    # first row     first row       last row       last row
    rowspec[1,1] = rowspec[1,2] = rowspec[n,3] = rowspec[n,4] = 1
    colspec = spzeros(T,4,n)
    # third col     fourth col     third-to-last col  fourth-to-last col
    colspec[1,3] = colspec[2,4] = colspec[3,n-2] = colspec[4,n-3] = 1
    valspec = zeros(T,4,4)
    # [1,3]          [n,n-2]
    valspec[1,1] = valspec[3,3] = 3
    # [1,4]          [n,n-3]
    valspec[2,2] = valspec[4,4] = -1

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), rowspec, valspec, colspec), zeros(TCoefs, n)
end

function prefiltering_system{T,TCoefs,GT<:GridType}(::Type{T}, ::Type{TCoefs}, n::Int, ::Type{Quadratic{Periodic}}, ::Type{GT})
    dl,d,du = inner_system_diags(T,n,Quadratic{Periodic})

    rowspec = spzeros(T,n,2)
    # first row       last row
    rowspec[1,1] = rowspec[n,2] = 1
    colspec = spzeros(T,2,n)
    # last col         first col
    colspec[1,n] = colspec[2,1] = 1
    valspec = zeros(T,2,2)
    # [1,n]            [n,1]
    valspec[1,1] = valspec[2,2] = SimpleRatio(1,8)

    Woodbury(lufact!(Tridiagonal(dl, d, du), Val{false}), rowspec, valspec, colspec), zeros(TCoefs, n)
end

sqr(x) = x*x
