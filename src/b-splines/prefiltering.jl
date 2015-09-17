deval{N}(::Val{N}) = N
padding{IT<:BSpline}(::Type{IT}) = Val{0}()
@generated function padding{IT}(t::Type{IT})
    pad = [deval(padding(IT.parameters[d])) for d = 1:length(IT.parameters)]
    t = tuple(pad...)
    :(Val{$t}())
end

function padded_index{N,pad}(sz::NTuple{N,Int}, ::Val{pad})
    szpad = ntuple(i->sz[i]+2padextract(pad,i), N)::NTuple{N,Int}
    ind = Array(UnitRange{Int},N)
    for i in 1:N
        p = padextract(pad,i)
        ind[i] = 1+p:szpad[i]-p
    end
    ind,szpad
end

copy_with_padding{IT}(A, ::Type{IT}) = copy_with_padding(eltype(A), A, IT)
function copy_with_padding{TCoefs,IT<:DimSpec{InterpolationType}}(::Type{TCoefs}, A, ::Type{IT})
    Pad = padding(IT)
    ind,sz = padded_index(size(A), Pad)
    if sz == size(A)
        coefs = copy!(similar(A,TCoefs), A)
    else
        coefs = zeros(TCoefs, sz...)
        coefs[ind...] = A
    end
    coefs, Pad
end

prefilter!{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = A
prefilter{TWeights, TCoefs, IT<:BSpline, GT<:GridType}(::Type{TWeights}, ::Type{TCoefs}, A, ::Type{IT}, ::Type{GT}) = prefilter!(TWeights, copy!(similar(A,TCoefs), A), IT, GT), Val{0}()

function prefilter{TWeights,TCoefs,TSrc,N,IT<:Quadratic,GT<:GridType}(
    ::Type{TWeights}, ::Type{TCoefs}, A::Array{TSrc,N}, ::Type{BSpline{IT}}, ::Type{GT}
    )
    ret, Pad = copy_with_padding(TCoefs,A, BSpline{IT})
    prefilter!(TWeights, ret, BSpline{IT}, GT), Pad
end

function prefilter{TWeights,TCoefs,TSrc,N,IT<:Tuple{Vararg{Union(BSpline,NoInterp)}},GT<:DimSpec{GridType}}(
    ::Type{TWeights}, ::Type{TCoefs}, A::Array{TSrc,N}, ::Type{IT}, ::Type{GT}
    )
    ret, Pad = copy_with_padding(TCoefs,A, IT)
    prefilter!(TWeights, ret, IT, GT), Pad
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

function prefilter!{TWeights,TCoefs,N,IT<:Tuple{Vararg{Union(BSpline,NoInterp)}},GT<:DimSpec{GridType}}(
    ::Type{TWeights}, ret::Array{TCoefs,N}, ::Type{IT}, ::Type{GT}
    )
    local buf, shape, retrs
    sz = size(ret)
    first = true
    for dim in 1:N
        it = iextract(IT, dim)
        if it != NoInterp
            M, b = prefiltering_system(TWeights, TCoefs, sz[dim], bsplinetype(it), iextract(GT, dim))
            if M != nothing
                A_ldiv_B_md!(ret, M, ret, dim, b)
            end
        end
    end
    ret
end

prefiltering_system(::Any, ::Any, ::Any, ::Any, ::Any) = nothing, nothing
