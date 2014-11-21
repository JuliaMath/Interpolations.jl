# Extrapolation behavior

## 1. Simple extrapolation

Extrapolation in `Interpolations.jl` takes several forms. The simplest, is when the extrapolation behavior can be described with the following pseudo-code:

    if (the coordinate is outside the domain) then
        do something that is
            * well defined without looking at the data set
            * decides the outcome (error/return value) of the indexing operation
    end

    proceed to inbounds interpolation

An example of this interpolation behavior is `ExtrapError`, which simply throws a bounds error if the coordinate is outside the domain. `Interpolations.jl` could support, for example, the following variants:

* `ExtrapError`: Throws a `BoundsError`, just like `Grid.jl`
* `ExtrapNaN`: Returns `convert(T, NaN)`, where `T` is `eltype(data)`
* `ExtrapNull`: Returns a value-less `Nullabel{T}`

## 2. Index transformation extrapolation

The next form is index transformation, which can be described as

    if (the coordinate is outside the domain) then
        calculate an index that is inside the domain, which gives
        the extrapolated value
    end

    proceed to inbounds interpolation (using the transformed index)

An example here is `ExtrapPeriodic`, which transforms the coordinate index to one which is inside the domain by means of modulo calculations. Another example is `ExtrapConstant`, which clamps out-of-bounds coordinates to their nearest inbounds data point.

For some of these, extra care needs to be taken in higher dimensions, when deciding what happens in the "outside corners":

              |               |     what happens here?
              |               |
    ----------+---------------+--------
              |               |
              |   the domain  |
              |               |
    ----------+---------------+--------
              |               |
              |               |    ...and here? 
