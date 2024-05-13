## ChainRulesCore Integration

The `WeightedIndex` infrastructure has issues working with autodiff libraries. Autodiff is facilitate by integration with ChainRulesCore. A custom [rrule](https://juliadiff.org/ChainRulesCore.jl/dev/index.html) is defined such that

```julia
y, itp_pullback = rrule(itp, 1)
```
`itp_pullback` takes a perturbation on `y` and returns how it effects each `x` dimension. Since `Interpolations` already has a `gradient` function, `pullback` reuses it by scaling it by `Δy`.

This enables integration with autodiff libraries like Zygote, enabling

```julia
x = 1:10
y = sin.(x)
itp = interpolate(y, BSpline(Linear()))
Zygote.gradient(itp, 2)
#([-0.7681774187658145],)
```
