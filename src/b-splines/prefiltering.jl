deval{N}(::Val{N}) = N
padding{IT<:BSpline}(::Type{IT}) = Val{0}()
@generated function padding{IT}(t::Type{IT})
    pad = [deval(padding(IT.parameters[d])) for d = 1:length(IT.parameters)]
    t = tuple(pad...)
    :(Val{$t}())
end

function padded_index{N,pad}(sz::NTuple{N,Int}, ::Val{pad})
    szpad = ntuple(i->sz[i]+2padextract(pad,i), N)::NTuple{N,Int}
    ind = Array{UnitRange{Int}}(N)
    for i in 1:N
        p = padextract(pad,i)
        ind[i] = 1+p:szpad[i]-p
    end
    ind,szpad
end

copy_with_padding{IT}(A, ::Type{IT}) = copy_with_padding(eltype(A), A, IT)
function copy_with_padding{TC,IT<:DimSpec{InterpolationType}}(::Type{TC}, A, ::Type{IT})
    Pad = padding(IT)
    ind,sz = padded_index(size(A), Pad)
    if sz == size(A)
        coefs = copy!(Array{TC}(size(A)), A)
    else
        coefs = zeros(TC, sz...)
        coefs[ind...] = A
    end
    coefs, Pad
end

prefilter!{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = A
prefilter{TWeights, TC, IT<:BSpline, GT<:GridType}(::Type{TWeights}, ::Type{TC}, A, ::Type{IT}, ::Type{GT}) = prefilter!(TWeights, copy!(Array{TC}(size(A)), A), IT, GT), Val{0}()

function prefilter{TWeights,TC,IT<:Union{Cubic,Quadratic},GT<:GridType}(
    ::Type{TWeights}, ::Type{TC}, A::AbstractArray, ::Type{BSpline{IT}}, ::Type{GT}
    )
    ret, Pad = copy_with_padding(TC, A, BSpline{IT})
    prefilter!(TWeights, ret, BSpline{IT}, GT), Pad
end

function prefilter{TWeights,TC,IT<:Tuple{Vararg{Union{BSpline,NoInterp}}},GT<:DimSpec{GridType}}(
    ::Type{TWeights}, ::Type{TC}, A::AbstractArray, ::Type{IT}, ::Type{GT}
    )
    ret, Pad = copy_with_padding(TC, A, IT)
    prefilter!(TWeights, ret, IT, GT), Pad
end

function prefilter!{TWeights,TCoefs<:AbstractArray,IT<:Union{Quadratic,Cubic},GT<:GridType}(
    ::Type{TWeights}, ret::TCoefs, ::Type{BSpline{IT}}, ::Type{GT}
    )
    local buf, shape, retrs
    sz = size(ret)
    first = true
    for dim in 1:ndims(ret)
        M, b = prefiltering_system(TWeights, eltype(TCoefs), sz[dim], IT, GT)
        A_ldiv_B_md!(ret, M, ret, dim, b)
    end
    ret
end

function prefilter!{TWeights,TCoefs<:AbstractArray,IT<:Tuple{Vararg{Union{BSpline,NoInterp}}},GT<:DimSpec{GridType}}(
    ::Type{TWeights}, ret::TCoefs, ::Type{IT}, ::Type{GT}
    )
    local buf, shape, retrs
    sz = size(ret)
    first = true
    for dim in 1:ndims(ret)
        it = iextract(IT, dim)
        if it != NoInterp
            M, b = prefiltering_system(TWeights, eltype(TCoefs), sz[dim], bsplinetype(it), iextract(GT, dim))
            if M != nothing
                A_ldiv_B_md!(ret, M, ret, dim, b)
            end
        end
    end
    ret
end

prefiltering_system(::Any, ::Any, ::Any, ::Any, ::Any) = nothing, nothing

"""
    M, b = prefiltering_system{T,TC,GT<:GridType,D<:Degree}m(::T, ::Type{TC}, n::Int, ::Type{D}, ::Type{GT})

Given element types (`T`, `TC`) and interpolation scheme
(`GT`, `D`) as well the number of rows in the data input (`n`), compute the
system used to prefilter spline coefficients. Boundary conditions determine the
values on the first and last rows.

Some of these boundary conditions require that these rows have off-tridiagonal
elements (e.g the `[1,3]` element of the matrix). To maintain the efficiency of
solving tridiagonal systems, the [Woodbury matrix identity](https://en.wikipedia.org/wiki/Woodbury_matrix_identity)
is used to add additional elements off the main 3 diagonals.

The filtered coefficients are given by solving the equation system

    M * c = v + b

where `c` are the sought coefficients, and `v` are the data points.
"""
prefiltering_system

"""
    dl, d, du = inner_system_diags{T,IT}(::Type{T}, n::Int, ::Type{IT})

Helper function to generate the prefiltering equation system: generates the diagonals
for a `n`-by-`n` tridiagonal matrix with eltype `T` corresponding to the interpolation
type `IT`.

`dl`, `d`, and `du` are intended to be used e.g. as in `M = Tridiagonal(dl, d, du)`
"""
inner_system_diags
