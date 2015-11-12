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
function copy_with_padding{TC,IT<:DimSpec{InterpolationType}}(::Type{TC}, A, ::Type{IT})
    Pad = padding(IT)
    ind,sz = padded_index(size(A), Pad)
    if sz == size(A)
        coefs = copy!(Array(TC, size(A)), A)
    else
        coefs = zeros(TC, sz...)
        coefs[ind...] = A
    end
    coefs, Pad
end

prefilter!{TWeights, IT<:BSpline, GT<:GridType}(::Type{TWeights}, A, ::Type{IT}, ::Type{GT}) = A
prefilter{TWeights, TC, IT<:BSpline, GT<:GridType}(::Type{TWeights}, ::Type{TC}, A, ::Type{IT}, ::Type{GT}) = prefilter!(TWeights, copy!(Array(TC, size(A)), A), IT, GT), Val{0}()

function prefilter{TWeights,TC,IT<:Quadratic,GT<:GridType}(
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

function prefilter!{TWeights,TCoefs<:AbstractArray,IT<:Quadratic,GT<:GridType}(
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
