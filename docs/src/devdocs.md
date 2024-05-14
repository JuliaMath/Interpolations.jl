```@meta
DocTestSetup = quote
    using Interpolations
end
```

# Developer documentation

Conceptually, `Interpolations.jl` supports two operations: construction and
usage of interpolants.

## Interpolant construction

Construction creates the interpolant object. In some situations this is
relatively trivial: for example, when using only `NoInterp`, `Constant`, or
`Linear` interpolation schemes, construction essentially corresponds to
recording the array of values and the "settings" (the interpolation scheme)
specified at the time of construction. This case is simple because interpolated
values may be efficiently computed directly from the on-grid values supplied at
construction time: ``(1-\Delta x) a_i  + \Delta x a_{i+1}`` reconstructs ``a_i``
when ``\Delta x = 0``.

For `Quadratic` and higher orders, efficient computation requires that the array
of values be *prefiltered*.
This essentially corresponds to "inverting" the computation that will be
performed during interpolation, so as to approximately reconstruct the original
values at on-grid points. Generally speaking this corresponds to solving a
nearly-tridiagonal system of equations, inverting an underlying interpolation
scheme such as
``p(\Delta x) \tilde a_{i-1} + q(\Delta x) \tilde a_i + p(1-\Delta x) \tilde a_{i+1}``
for some functions ``p`` and ``q`` (see [`Quadratic`](@ref) for further details).
Here ``\tilde a`` is the pre-filtered version of `a`, designed so that
substituting ``\Delta x = 0`` (for which one may *not* get 0 and 1 for the
``p`` and ``q`` calls, respectively) approximately recapitulates ``a_i``.

The exact system of equations to be solved depends on the interpolation order
and boundary conditions.
Boundary conditions often introduce deviations from perfect tridiagonality;
these "extras" are handled efficiently by the `WoodburyMatrices` package.
These computations are implemented independently along each axis using the
`AxisAlgorithms` package.

In the `doc` directory there are some old files that give some of the
mathematical details.  A useful reference is [Thévenaz2000](@citet).

!!! note
    As an application of these concepts, note that supporting quadratic or cubic
    interpolation for `Gridded` would only require that someone implement
    prefiltering schemes for non-uniform grids;
    it's just a question of working out a little bit of math.

## Interpolant usage

Usage occurs when evaluating `itp(x, y...)`, or `Interpolations.gradient(itp, x, y...)`,
etc.  Usage itself involves two sub-steps: computation of the *weights* and then
performing the *interpolation*.

### Weight computation

Weights depend on the interpolation scheme and the location `x, y...` but
*not* the coefficients of the array we are interpolating.
Consequently there are many circumstances where one might want to reuse
previously-computed weights, and `Interpolations.jl` has been carefully
designed with that kind of reuse in mind.

The key concept here is the [`Interpolations.WeightedIndex`](@ref), and there is
no point repeating its detailed docstring here.
It suffices to add that `WeightedIndex` is actually an abstract type, with two
concrete subtypes:

- `WeightedAdjIndex` is for indexes that will address *adjacent* points of the
  coefficient array (ones where the index increments by 1 along the corresponding
  dimension). These are used when prefiltering produces padding that can be used
  even at the edges, or for schemes like `Linear` interpolation which require no
  padding.
- `WeightedArbIndex` stores both the weight and index associated with each
  accessed grid point, and can therefore encode grid access patterns. These are
  used in specific circumstances--a prime example being periodic boundary
  conditions--where the coefficients array may be accessed at something other than
  adjacent locations.

`WeightedIndex` computation reflects the interpolation scheme (e.g., `Linear` or `Quadratic`) and also whether one is computing values, `gradient`s, or `hessian`s. The handling of derivatives will be described further below.

### Interpolation

General `AbstractArray`s may be indexed with `WeightedIndex` indices,
and the result produces the interpolated value. In other words, the end result
is just `itp.coefs[wis...]`, where `wis` is a tuple of `WeightedIndex` indices.
To make sure this overloading is effective, we wrap the `coefs` with `InterpGetindex`，
i.e. `InterpGetindex(itp.coefs)[wis...]`

Derivatives along a particular axis can be computed just by substituting a component of `wis` for one that has been designed to compute derivatives rather than values.

As a demonstration, let's see how the following computation occurs:

```jldoctest derivs
julia> A = reshape(1:27, 3, 3, 3)
3×3×3 reshape(::UnitRange{Int64}, 3, 3, 3) with eltype Int64:
[:, :, 1] =
 1  4  7
 2  5  8
 3  6  9

[:, :, 2] =
 10  13  16
 11  14  17
 12  15  18

[:, :, 3] =
 19  22  25
 20  23  26
 21  24  27

julia> itp = interpolate(A, BSpline(Linear()));

julia> x = (1.2, 1.4, 1.7)
(1.2, 1.4, 1.7)

julia> itp(x...)
8.7
```

!!! note
    By using the debugging facilities of an IDE like Juno or VSCode,
    or using Debugger.jl from the REPL, you can easily step in to
    the call above and follow along with the description below.

Aside from details such as bounds-checking, the key call is to
[`Interpolations.weightedindexes`](@ref):

```jldoctest derivs
julia> wis = Interpolations.weightedindexes((Interpolations.value_weights,), Interpolations.itpinfo(itp)..., x)
(Interpolations.WeightedAdjIndex{2, Float64}(1, (0.8, 0.19999999999999996)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.6000000000000001, 0.3999999999999999)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.30000000000000004, 0.7)))

julia> wis[1]
Interpolations.WeightedAdjIndex{2, Float64}(1, (0.8, 0.19999999999999996))

julia> wis[2]
Interpolations.WeightedAdjIndex{2, Float64}(1, (0.6000000000000001, 0.3999999999999999))

julia> wis[3]
Interpolations.WeightedAdjIndex{2, Float64}(1, (0.30000000000000004, 0.7))

julia> Interpolations.InterpGetindex(A)[wis...]
8.7
```

You can see that each of `wis` corresponds to a specific position: 1.2, 1.4, and 1.7 respectively. We can index `A` at `wis`, and it returns the value of `itp(x...)`, which here is just
```plain
  0.8 * A[1, wis[2], wis[3]] + 0.2 * A[2, wis[2], wis[3]]
= 0.6 * (0.8 * A[1, 1, wis[3]] + 0.2 * A[2, 1, wis[3]]) +
  0.4 * (0.8 * A[1, 2, wis[3]] + 0.2 * A[2, 2, wis[3]])
= 0.3 * (0.6 * (0.8 * A[1, 1, 1] + 0.2 * A[2, 1, 1]) +
         0.4 * (0.8 * A[1, 2, 1] + 0.2 * A[2, 2, 1])  ) +
  0.7 * (0.6 * (0.8 * A[1, 1, 2] + 0.2 * A[2, 1, 2]) +
         0.4 * (0.8 * A[1, 2, 2] + 0.2 * A[2, 2, 2])  )
```
This computed the value of `itp` at `x...` because we called `weightedindexes`
with just a single function, [`Interpolations.value_weights`](@ref) (meaning,
"the weights needed to compute the value").

!!! note
    Remember that prefiltering is not used for `Linear` interpolation.
    In a case where prefiltering is used, we would substitute
    `InterpGetindex(itp.coefs)[wis...]` for `InterpGetindex(A)[wis...]` above.

To compute derivatives, we *also* pass additional functions like [`Interpolations.gradient_weights`](@ref):

```jldoctest derivs
julia> wis = Interpolations.weightedindexes((Interpolations.value_weights, Interpolations.gradient_weights), Interpolations.itpinfo(itp)..., x)
((Interpolations.WeightedAdjIndex{2, Float64}(1, (-1.0, 1.0)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.6000000000000001, 0.3999999999999999)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.30000000000000004, 0.7))), (Interpolations.WeightedAdjIndex{2, Float64}(1, (0.8, 0.19999999999999996)), Interpolations.WeightedAdjIndex{2, Float64}(1, (-1.0, 1.0)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.30000000000000004, 0.7))), (Interpolations.WeightedAdjIndex{2, Float64}(1, (0.8, 0.19999999999999996)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.6000000000000001, 0.3999999999999999)), Interpolations.WeightedAdjIndex{2, Float64}(1, (-1.0, 1.0))))

julia> wis[1]
(Interpolations.WeightedAdjIndex{2, Float64}(1, (-1.0, 1.0)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.6000000000000001, 0.3999999999999999)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.30000000000000004, 0.7)))

julia> wis[2]
(Interpolations.WeightedAdjIndex{2, Float64}(1, (0.8, 0.19999999999999996)), Interpolations.WeightedAdjIndex{2, Float64}(1, (-1.0, 1.0)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.30000000000000004, 0.7)))

julia> wis[3]
(Interpolations.WeightedAdjIndex{2, Float64}(1, (0.8, 0.19999999999999996)), Interpolations.WeightedAdjIndex{2, Float64}(1, (0.6000000000000001, 0.3999999999999999)), Interpolations.WeightedAdjIndex{2, Float64}(1, (-1.0, 1.0)))

julia> Interpolations.InterpGetindex(A)[wis[1]...]
1.0

julia> Interpolations.InterpGetindex(A)[wis[2]...]
3.000000000000001

julia> Interpolations.InterpGetindex(A)[wis[3]...]
9.0
```
In this case you can see that `wis` is a 3-tuple-of-3-tuples. `A[wis[i]...]` can
be used to compute the `i`th component of the gradient.

If you look carefully at each of the entries in `wis`, you'll see that each
"inner" 3-tuple copies two of the three elements in the `wis` we obtained when
we called `weightedindexes` with just `value_weights` above.  `wis[1]` replaces
the first entry with a weighted index having weights `(-1.0, 1.0)`, which
corresponds to computing the *slope* along this dimension.
Likewise `wis[2]` and `wis[3]` replace the second and third value-index,
respectively, with the same slope computation.

Hessian computation is quite similar, with the difference that one sometimes
needs to replace two different indices or the same index with a set of weights
corresponding to a second derivative.

Consequently derivatives along particular directions are computed simply by
"weight replacement" along the corresponding dimensions.

The code to do this replacement is a bit complicated due to the need to support
arbitrary dimensionality in a manner that allows Julia's type-inference to succeed.
It makes good use of *tuple* manipulations, sometimes called "lispy tuple programming."
You can search Julia's discourse forum for more tips about how to program this way.
It could alternatively be done using generated functions, but this would
increase compile time considerably and can lead to world-age problems.

### GPU Support
At present, `Interpolations.jl` supports interpolant usage on GPU via broadcasting.

A basic work flow looks like:
```julia
using Interpolations, Adapt, CUDA # Or any other GPU package
itp = Interpolations.interpolate([1, 2, 3], BSpline(Linear())); # construct the interpolant object on CPU
cuitp = adapt(CuArray{Float32}, itp); # adapt it to GPU memory
cuitp.(1:0.5:2) # call interpolant object via broadcast
gradient.(Ref(cuitp), 1:0.5:2)
```

To achieve this, an `ITP <: AbstractInterpolation` should define it's own
`Adapt.adapt_structure(to, itp::ITP)`, which constructs a new `ITP` with the
adapted fields (`adapt(to, itp.fieldname)`) of `itp`. The field adaption could
be skipped if we know that it has been GPU-compatable, e.g. a `isbit` range.

!!! note
    Some adaptors may change the storage type. Please ensure that the adapted
    `itp` has the correct element type via the method `eltype`.

Also, all GPU-compatable `AbstractInterpolation`s should define their own
`Interpolations.root_storage_type`. This function allows us to modify the
broadcast mechanism by overloading the default `BroadcastStyle`. See
[Customizing broadcasting](https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting)
for more details.


```@bibliography
Pages = ["devdocs.md"]
Canonical = false
```
