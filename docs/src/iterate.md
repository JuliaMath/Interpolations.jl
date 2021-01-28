# Knot Iteration

Given an `AbstractInterpolation` `itp` get an iterator over it's knots using
`knots(itp)`

```julia
using Interpolations
itp = interpolate(rand(4), options...)
kiter = knots(itp) # Iterator over knots
collect(kiter) # Array of knots [1, 2, 3, 4]

```

For multiple dimensions, the iterator will return tuples of positions
(ie. `(x, y, ...)`),  with the first coordinate changing the fastest.

```jldoctest iterate-interpolate; setup = :(using Interpolations)
julia> itp = interpolate(ones(3,3), BSpline(Linear()));

julia> kiter = knots(itp);

julia> collect(kiter)
3×3 Array{Tuple{Int64,Int64},2}:
 (1, 1)  (1, 2)  (1, 3)
 (2, 1)  (2, 2)  (2, 3)
 (3, 1)  (3, 2)  (3, 3)
```

The number of elements and size of the iterator can be found as shown:

```jldoctest iterate-interpolate; setup = :(using Interpolations)
julia> length(kiter)
9

julia> size(kiter)
(3, 3)

```


## Extrapolated Knots

Given an `AbstractExtrapolation` `etp`, `knots(etp)` will also iterate over the
the knots with the following behavior.

- For `Throw`, `Flat`, `Line` iterate the knots once
- For `Periodic` and `Reflect` generate an infinite sequence of knots starting
  at the first knot.

As `Periodic` and `Reflect` generate infinite sequences of knots, `length` and
`size` are undefined. For `Throw`, `Flat`, `Line`, `length` and `size` behave as
expected.

### Periodic

With Periodic boundary condition, knots repeat indefinitely with the first and
last knot being co-located. (ie. in the example below `etp(2.0) = 1.0` not
`4.0`).

```jldoctest periodic-demo; setup = :(using Interpolations)
julia> x = [1.0, 1.5, 1.75, 2.0];

julia> etp = LinearInterpolation(x, x.^2, extrapolation_bc=Periodic());

julia> kiter = knots(etp);

julia> k = Iterators.take(kiter, 6) |> collect
6-element Array{Float64,1}:
 1.0
 1.5
 1.75
 2.0
 2.5
 2.75

```

Extrapolating to the generated knots `etp.(k)`, confirms that the extrapolated
knots do map back to the correct inbound knots (ie. `etp(k[1]) == etp(k[4])`).

```jldoctest periodic-demo
julia> etp.(k)
6-element Array{Float64,1}:
 1.0
 2.25
 3.0625
 1.0
 2.25
 3.0625

```

### Reflect

With the `Reflect` boundary condition knots repeat indefinitely, following the
pattern shown below (Offset terms are not shown for brevity).

```
k[1], k[2], ..., k[end-1], k[end], k[end+1], ... k[2], k[1], k[2], ...
```

```jldoctest reflect-demo; setup = :(using Interpolations)
julia> x = [1.0, 1.5, 1.75, 2.0];

julia> etp = LinearInterpolation(x, x.^2, extrapolation_bc=Reflect());

julia> kiter = knots(etp);

julia> k = Iterators.take(kiter, 6) |> collect
6-element Array{Float64,1}:
 1.0
 1.5
 1.75
 2.0
 2.25
 2.5

```

Evaluating the extrapolation at `etp.(k)` confirms that the extrapolated knots
correspond to the correct inbound knots.

```jldoctest reflect-demo
julia> etp.(k)
6-element Array{Float64,1}:
 1.0
 2.25
 3.0625
 4.0
 3.0625
 2.25

```

### Multiple Dimensions

As with an `AbstractInterpolation`, iterating over knots for a
multi-dimensional extrapolation also supported.

```jldoctest; setup = :(using Interpolations)
julia> x = [1.0, 1.5, 1.75, 2.0];

julia> etp = LinearInterpolation((x, x), x.*x');

julia> knots(etp) |> collect
4×4 Array{Tuple{Float64,Float64},2}:
 (1.0, 1.0)   (1.0, 1.5)   (1.0, 1.75)   (1.0, 2.0)
 (1.5, 1.0)   (1.5, 1.5)   (1.5, 1.75)   (1.5, 2.0)
 (1.75, 1.0)  (1.75, 1.5)  (1.75, 1.75)  (1.75, 2.0)
 (2.0, 1.0)   (2.0, 1.5)   (2.0, 1.75)   (2.0, 2.0)

```

Because some boundary conditions generate an infinite sequences of knots,
iteration over knots can end up "stuck" iterating along a single axis:

```jldoctest; setup = :(using Interpolations)
julia> x = [1.0, 1.5, 1.75, 2.0];

julia> etp = LinearInterpolation((x, x), x.*x', extrapolation_bc=(Periodic(), Throw()));

julia> knots(etp) |> k -> Iterators.take(k, 6) |> collect
6-element Array{Tuple{Float64,Float64},1}:
 (1.0, 1.0)
 (1.5, 1.0)
 (1.75, 1.0)
 (2.0, 1.0)
 (2.5, 1.0)
 (2.75, 1.0)

```

Rearranging the axes so non-repeating knots are first can address this issue:

```jldoctest; setup = :(using Interpolations)
julia> x = [1.0, 1.5, 1.75, 2.0];

julia> etp = LinearInterpolation((x, x), x.*x', extrapolation_bc=(Throw(), Periodic()));

julia> knots(etp) |> k -> Iterators.take(k, 6) |> collect
6-element Array{Tuple{Float64,Float64},1}:
 (1.0, 1.0)
 (1.5, 1.0)
 (1.75, 1.0)
 (2.0, 1.0)
 (1.0, 1.5)
 (1.5, 1.5)

```

### Directional Boundary Conditions

If the boundary conditions are directional, the forward boundary condition is
used to determine if the iterator will generate an infinite sequence of knots.

For example the following extrapolation `etp`, will throw an error for values
less than `1.0`, but will use `Periodic` extrapolation for values above `2.0`. As a
result, the iterator will generate an infinite sequence of knots starting at `1.0`.

```jldoctest iterate-directional-unbounded; setup = :(using Interpolations)
julia> x = [1.0, 1.2, 1.3, 2.0];

julia> etp = LinearInterpolation(x, x.^2, extrapolation_bc=((Throw(), Periodic()),));

julia> kiter = knots(etp);

julia> kiter |> k -> Iterators.take(k, 5) |> collect
5-element Array{Float64,1}:
 1.0
 1.2
 1.3
 2.0
 2.2

```

We can also check if the iterator has a length using: `Base.IteratorSize`

```jldoctest iterate-directional-unbounded
julia> Base.IteratorSize(kiter)
Base.IsInfinite()

```

Swapping the boundary conditions, results in a finite sequence of knots from
`1.0` to `2.0`.

```jldoctest iterate-directional-bounded; setup = :(using Interpolations)
julia> x = [1.0, 1.2, 1.3, 2.0];

julia> etp = LinearInterpolation(x, x.^2, extrapolation_bc=((Periodic(), Throw()),));

julia> kiter = knots(etp);

julia> collect(kiter)
4-element Array{Float64,1}:
 1.0
 1.2
 1.3
 2.0

```

As expected the iterator now has a defined length:

```jldoctest iterate-directional-bounded; setup = :(using Interpolations)
julia> Base.IteratorSize(kiter)
Base.HasLength()

julia> length(kiter)
4

julia> size(kiter)
(4,)

```

## `knotsbetween`

Given an `AbstractInterpolation` `itp` or results of `knots(itp)` get an iterator
over knots between `start` and `stop`.

```julia
using Interpolations
itp = interpolate(rand(4), options...)
krange = knotsbetween(itp; start=1.2, stop=3.0)
collect(kiter) # Array of knots between 1.2 and 3.0

```

We can iterate over all knots greater than `start` by omitting `stop`

```jldoctest knotsbetween-usage; setup = :(using Interpolations)
julia> x = [1.0, 1.5, 1.75, 2.0];

julia> etp = LinearInterpolation(x, x.^2, extrapolation_bc=Periodic());

julia> krange = knotsbetween(etp; start=4.0);

julia> Iterators.take(krange, 5) |> collect
5-element Array{Float64,1}:
 4.5
 4.75
 5.0
 5.5
 5.75

```

If we omit `start`, iteration will range from the first knot to `stop`

```jldoctest knotsbetween-usage
julia> krange = knotsbetween(etp; stop=4.0);

julia> collect(krange)
9-element Array{Float64,1}:
 1.0
 1.5
 1.75
 2.0
 2.5
 2.75
 3.0
 3.5
 3.75

```

It is an error to not provided `start` and `stop`

```jldoctest knotsbetween-usage
julia> knotsbetween(etp)
ERROR: ArgumentError: At least one of `start` or `stop` must be specified
[...]
```

### Multiple Dimensions

When used with a multi-dimensional interpolant, `knotsbetween` can be used to
iterate overall knots such that: `start[i] < k[i] stop[i]` where `i` indexes
dimensions.

```jldoctest; setup = :(using Interpolations)
julia> x = [1.0, 1.5, 1.75, 2.0];

julia> etp = LinearInterpolation((x, x), x.*x');

julia> knotsbetween(etp; start=(1.2, 1.5), stop=(1.8, 3.0)) |> collect
2×2 Array{Tuple{Float64,Float64},2}:
 (1.5, 1.75)   (1.5, 2.0)
 (1.75, 1.75)  (1.75, 2.0)

```
