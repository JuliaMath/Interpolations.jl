padding{IT<:BSpline}(::Type{IT}) = 0

function padded_index{N}(sz::NTuple{N,Int}, pad)
    sz = [s+2pad for s in sz]
    ind = Array(Any,N)
    for i in 1:N
        ind[i] = 1+pad:sz[i]-pad
    end
    ind,sz
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

@generated function prefilter{TWeights,TCoefs,N,IT<:Quadratic,GT<:GridType}(
    ::Type{TWeights}, A::Array{TCoefs,N}, ::Type{BSpline{IT}}, ::Type{GT}
    )
    quote
        ret, pad = copy_with_padding(A, BSpline{IT})

        szs = collect(size(ret))
        strds = collect(strides(ret))

        for dim in 1:$N
            n = szs[dim]
            szs[dim] = 1

            M, b = prefiltering_system(TWeights, TCoefs, n, IT, GT)

            @nloops $N i d->1:szs[d] begin
                cc = @ntuple $N i

                strt_diff = sum([(cc[i]-1)*strds[i] for i in 1:length(cc)])
                strt = 1 + strt_diff
                rng = range(strt, strds[dim], n)

                bdiff = ret[rng]
                b += bdiff
                ret[rng] = M \ b
                b -= bdiff
            end
            szs[dim] = n
        end
        ret
    end
end
