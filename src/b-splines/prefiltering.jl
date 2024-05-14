padded_axes(axs, it::InterpolationType) = (ax->padded_axis(ax, it)).(axs)
padded_axes(axs::NTuple{N,AbstractUnitRange}, it::NTuple{N,InterpolationType}) where N =
    padded_axis.(axs, it)

padded_similar(::Type{TC}, inds::Tuple{Vararg{Base.OneTo{Int}}}) where TC = Array{TC}(undef, length.(inds))
padded_similar(::Type{TC}, inds) where TC = OffsetArray{TC}(undef, inds)

# Narrow ax by the amount that axpad is larger
padinset(ax::AbstractUnitRange, axpad) = 2*first(ax)-first(axpad):2*last(ax)-last(axpad)
function padinset(ax::Base.OneTo, axpad::Base.OneTo)
    # We don't have any types that pad asymmetrically. Therefore if they both start at 1,
    # they must be the same
    @assert ax == axpad
    return ax
end

ct!(coefs, indscp, A, indsA) = copyto!(coefs, CartesianIndices(indscp), A, CartesianIndices(indsA))

copy_with_padding(A, it) = copy_with_padding(eltype(A), A, it)
function copy_with_padding(::Type{TC}, A, it::DimSpec{InterpolationType}) where {TC}
    indsA = axes(A)
    indspad = padded_axes(indsA, it)
    coefs = padded_similar(TC, indspad)
    if indspad == indsA
        coefs = copyto!(coefs, A)
    else
        fill!(coefs, zero(TC))
        ct!(coefs, indsA, A, indsA)
    end
    coefs
end

prefilter!(::Type{TWeights}, A::AbstractArray, ::BSpline{D}, ::GridType) where {TWeights,D<:Union{Constant,Linear}} = A

function prefilter(
    ::Type{TWeights}, ::Type{TC}, A::AbstractArray,
    it::Union{BSpline,Tuple{Vararg{Union{BSpline,NoInterp}}}}
    ) where {TWeights,TC}
    ret = copy_with_padding(TC, A, it)
    prefilter!(TWeights, ret, it)
end

function prefilter(
    ::Type{TWeights}, ::Type{TC}, A::AbstractArray,
    it::Union{BSpline,Tuple{Vararg{Union{BSpline,NoInterp}}}},
    λ::Real, k::Int
    ) where {TWeights,TC}
    ret = copy_with_padding(TC, A, it)
    prefilter!(TWeights, ret, it, λ, k)
end

function prefilter!(
    ::Type{TWeights}, ret::TCoefs, it::BSpline
    ) where {TWeights,TCoefs<:AbstractArray}
    sz = size(ret)
    first = true
    for dim in 1:ndims(ret)
        M, b = prefiltering_system(TWeights, eltype(TCoefs), sz[dim], degree(it))
        A_ldiv_B_md!(popwrapper(ret), M, popwrapper(ret), dim, b)
    end
    ret
end

function diffop(::Type{TWeights}, n::Int, k::Int) where {TWeights}
    D = spdiagm(0 => ones(TWeights,n))
    for i in 1:k
        D = diff(D; dims=1)
    end
    ## TODO: Normalize by n?
    D' * D
end
diffop(n::Int, k::Int) = diffop(Float64, n, k)
### TODO: add compiled constructor for most common operators of order k=1,2


function prefilter!(
    ::Type{TWeights}, ret::TCoefs, it::BSpline, λ::Real, k::Int
    ) where {TWeights,TCoefs<:AbstractArray}

    sz = size(ret)

    # Test if regularizing is allowed
    fallback = if ndims(ret) > 1
        @warn "Smooth BSpline only available for Vectors, fallback to non-smoothed"
        true
    elseif λ < 0
        @warn "Smooth BSpline require a non-negative λ, fallback to non-smoothed: $(λ)"
        true
    elseif λ == 0
        true
    else
        false
    end

    if fallback
        return prefilter!(TWeights, ret, it)
    end

    M, b = prefiltering_system(TWeights, eltype(TCoefs), sz[1], degree(it))
    ## Solve with regularization
    # Convert to dense Matrix because `*(Woodbury{Adjoint}, Woodbury)` is not defined
    Mt = Matrix(M)'
    Q = diffop(TWeights, sz[1], k)
    K = cholesky(Mt * Matrix(M) + λ * Q)
    B = Mt * popwrapper(ret)
    ldiv!(popwrapper(ret), K, B)

    ret
end

function prefilter!(
    ::Type{TWeights}, ret::TCoefs, its::Tuple{Vararg{Union{BSpline,NoInterp}}}
    ) where {TWeights,TCoefs<:AbstractArray}
    sz = size(ret)
    for dim in 1:ndims(ret)
        it = iextract(its, dim)
        if it != NoInterp
            M, b = prefiltering_system(TWeights, eltype(TCoefs), sz[dim], degree(it))
            if M != nothing
                A_ldiv_B_md!(popwrapper(ret), M, popwrapper(ret), dim, b)
            end
        end
    end
    ret
end

prefiltering_system(::Any, ::Any, ::Any, ::Any) = nothing, nothing

popwrapper(A) = A
popwrapper(A::OffsetArray) = A.parent

"""
    M, b = prefiltering_system{T,TC,GT<:GridType,D<:Degree}(::T, ::Type{TC}, n::Int, ::Type{D}, ::Type{GT})

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
