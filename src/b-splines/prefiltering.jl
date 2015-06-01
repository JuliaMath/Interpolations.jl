padding{IT<:BSpline}(::Type{IT}) = Val{0}()

function padded_index{N,pad}(sz::NTuple{N,Int}, ::Val{pad})
    szpad = ntuple(i->sz[i]+2pad, N)::NTuple{N,Int}
    ind = Array(UnitRange{Int},N)
    for i in 1:N
        ind[i] = 1+pad:szpad[i]-pad
    end
    ind,szpad
end
function copy_with_padding{IT<:InterpolationType}(A, ::Type{IT})
    Pad = padding(IT)
    ind,sz = padded_index(size(A), Pad)
    coefs = zeros(eltype(A), sz...)
    coefs[ind...] = A
    coefs, Pad
end

prefilter!{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = A
prefilter{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = prefilter!(TWeights, copy(A), IT, GT), Val{0}()

function prefilter{TWeights,TCoefs,N,IT<:Quadratic,GT<:GridType}(
    ::Type{TWeights}, A::Array{TCoefs,N}, ::Type{BSpline{IT}}, ::Type{GT}
    )
    ret, Pad = copy_with_padding(A, BSpline{IT})
    prefilter!(TWeights, ret, BSpline{IT}, GT), Pad
end

function prefilter!{TWeights,TCoefs,N,IT<:Quadratic,GT<:GridType}(
    ::Type{TWeights}, ret::Array{TCoefs,N}, ::Type{BSpline{IT}}, ::Type{GT}
    )
    local buf, shape, retrs
    sz = size(ret)
    first = true
    for dim in 1:N
        M, b = prefiltering_system(TWeights, TCoefs, sz[dim], IT, GT)
        A_ldiv_B_md!(ret, M, ret, dim, b)
    end
    ret
end
