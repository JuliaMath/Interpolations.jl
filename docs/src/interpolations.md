## General usage

Note: the current version of `Interpolations` supports interpolation evaluation using index calls `[]`, but this feature will be deprecated in future. We highly recommend function calls with `()` as follows.

Given an `AbstractArray` `A`, construct an "interpolation object" `itp` as
```julia
itp = interpolate(A, options...)
```
where `options...` (discussed below) controls the type of
interpolation you want to perform.  This syntax assumes that the
samples in `A` are equally-spaced.

To evaluate the interpolation at position `(x, y, ...)`, simply do
```julia
v = itp(x, y, ...)
```

Some interpolation objects support computation of the gradient, which
can be obtained as
```julia
g = Interpolations.gradient(itp, x, y, ...)
```
or as
```julia
Interpolations.gradient!(g, itp, x, y, ...)
```
where `g` is a pre-allocated vector.

Some interpolation objects support computation of the hessian, which
can be obtained as
```julia
h = Interpolations.hessian(itp, x, y, ...)
```
or
```julia
Interpolations.hessian!(h, itp, x, y, ...)
```
where `h` is a pre-allocated matrix.

`A` may have any element type that supports the operations of addition
and multiplication.  Examples include scalars like `Float64`, `Int`,
and `Rational`, but also multi-valued types like `RGB` color vectors.

Positions `(x, y, ...)` are n-tuples of numbers. Typically these will
be real-valued (not necessarily integer-valued), but can also be of types
such as [DualNumbers](https://github.com/JuliaDiff/DualNumbers.jl) if
you want to verify the computed value of gradients.
(Alternatively, verify gradients using [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl).)
You can also use
Julia's iterator objects, e.g.,

```julia
function ongrid!(dest, itp)
    for I in CartesianIndices(itp)
        dest[I] = itp(I)
    end
end
```
would store the on-grid value at each grid point of `itp` in the output `dest`.
Finally, courtesy of Julia's indexing rules, you can also use
```julia
fine = itp(range(1,stop=10,length=1001), range(1,stop=15,length=201))
```

There is also an abbreviated [Convenience notation](@ref).
