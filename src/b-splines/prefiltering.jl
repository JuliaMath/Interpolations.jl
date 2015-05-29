padding{IT<:BSpline}(::Type{IT}) = 0

function padded_index{N}(sz::NTuple{N,Int}, pad)
    szpad = ntuple(i->sz[i]+2pad, N)::NTuple{N,Int}
    ind = Array(UnitRange{Int},N)
    for i in 1:N
        ind[i] = 1+pad:szpad[i]-pad
    end
    ind,szpad
end
function copy_with_padding{IT<:InterpolationType}(A, ::Type{IT})
    pad = padding(IT)
    ind,sz = padded_index(size(A), pad)
    coefs = zeros(eltype(A), sz...)
    coefs[ind...] = A
    coefs, pad
end

prefilter!{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = A
prefilter{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = prefilter!(TWeights, copy(A), IT, GT)

function prefilter{TWeights,TCoefs,N,IT<:Quadratic,GT<:GridType}(
    ::Type{TWeights}, A::Array{TCoefs,N}, ::Type{BSpline{IT}}, ::Type{GT}
    )
    ret, pad = copy_with_padding(A, BSpline{IT})
    sz = size(ret)
    first = true
    for dim in 1:N
        M, b = prefiltering_system(TWeights, TCoefs, sz[dim], IT, GT)
        if !isa(M, Woodbury)
            A_ldiv_B_md!(ret, M, ret, dim, b)
        else
            if first
                buf = Array(eltype(ret), length(ret,))
                shape = sz
                retrs = reshape(ret, shape)  # for type-stability against a future julia #10507
                first = false
            end
            bufrs = reshape(buf, shape)
            filter_dim1!(bufrs, M, retrs, b)
            shape = (sz[dim+1:end]..., sz[1:dim]...)::NTuple{N,Int}
            retrs = reshape(ret, shape)
            permutedims!(retrs, bufrs, ((2:N)..., 1))
        end
    end
    ret
end
