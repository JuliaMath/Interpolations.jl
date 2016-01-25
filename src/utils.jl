@inline sqr(x) = x*x
@inline cub(x) = x*x*x

"""
Build `rowspec`, `valspec`, `colspec` such that the product

`out = rowspec * valspec * colspec` will be equivalent to:

```julia
out = zeros(n, n)

for (i, j, v) in args
    out[i, j] = v
end
```

"""
function _build_woodbury_specs{T}(::Type{T}, n::Int, args::Tuple{Int, Int, Any}...)
    m = length(args)
    rowspec = spzeros(T, n, m)
    colspec = spzeros(T, m, n)
    valspec = zeros(T, m, m)

    ix = 1
    for (i, (row, col, val)) in enumerate(args)
        rowspec[row, ix] = 1
        colspec[ix, col] = 1
        valspec[ix, ix] = val
        ix += 1
    end

    rowspec, valspec, colspec
end
