# The math behind Interpolations.jl

Interpolations.jl uses [B-splines](http://en.wikipedia.org/wiki/B-spline#Definition), or *basis splines*, to interpolate a data set with a specified degree. In short, the idea is to use a set of basis functions with limited support, and in each interval between two data points, construct a linear combination of the basis functions which have support in that interval, each yielding a piece in a piecewise defined function with support across the entire data set. B-splines are defined to have the maximum degree of continuity, with minimal support. Specific details on the general mathematics of such interpolations can [this paper](http://dx.doi.org/10.1109/42.875199).

## The purpose of this document

There are a few core concepts one must understand in order to use `Interpolations.jl`:

1. Interpolation degree
2. Boundary conditions
3. Mid-point and on-grid interpolation
3. Extrapolation

Defining these three concepts, and how `Interpolations.jl` reasons about them, is the main focus of this document.

## 1. Interpolation degree

The interpolation degree decides what continuity properties the interpolant (i.e. the object which represents the interpolated data set) has; for example, a linear interpolation is a piecewise linear function, which is continuous but has discontinuous derivatives (or gradients, in higher dimensions).

Because of how B-splines are defined, the interpolation degree also decides the *support* of the *spline functions*. For a linear interpolation, the function value is decided by the value at *two* data points, so its support is two<sup>1.</sup>.

## 2. Boundary conditions

For higher interpolation degrees (specifically, from quadratic interpolation and up), the support of the B-splines is too large for the interpolating scheme to be able to figure out the values near the edges of the data sets. In order to close the equation systems at the edges, boundary conditions are used.

For quadratic interpolation, for example, a common boundary condition is to assume that the function is flat at the edges (i.e. the derivative there is 0). This lets us introduce an extra equation at each edge, through a finite approximation of the derivative, which closes the system. Another common way of terminating the interpolation is to extend the second-to-outermost all the way to the edge of the data set.

One subtlety concerns the location at which the boundary conditions are applied: at the edge grid point (`OnGrid()`) or at the halfway mark to the first beyond-the-edge index (`OnCell()`). `Interpolations.jl` supports both of these for interpolation schemes affected by boundary conditions (quadratic and cubic).

## Interlude: the `Interpolation` type hierarchy

All the above concepts fundamentally affect the behavior of an interpolating function, and in `Interpolations.jl` they therefore each have a representation in the type hierarchy.

Different types of interpolations are separated by *type parameters*, and each type of interpolation is a concrete type descending from

    abstract InterpolationType{D<:Degree,BC<:BoundaryCondition,G<:GridRepresentation}

For example, a quadratic on-grid implementation with flat boundary conditions, is represented by an `Interpolation{Quadratic,Flat,OnGrid}`, where `Quadratic`, `Flat` and `OnGrid` are in turn concrete types implementing the abstract types indicated above.

## 4. Extrapolation behavior

Somewhat orthogonal<sup>2</sup> to the concepts outlined above, is the concept of *extrapolation*, i.e. evaluation of the interpolant outside the domain defined by the data set. For some types of extrapolation, this behavior is defined by a translation of the interpolation coordinate to somewhere inside the domain (e.g. periodic or reflecting boundaries), while for other types it entails a separate calculation entirely (e.g. linear extrapolation).

[Read more](/doc/Extrapolation.md)

## Supporting an `Interpolation` type in `Interpolations.jl`

Not all combinations of the above four concepts are supported in the library; for example, neither constant nor linear on-grid interpolations need andy boundary conditions, and therefore they don't support any.

An interpolation is represented by an object of a concrete type descending from

    abstract Interpolation{IT<:InterpolationType,EB<:ExtrapolationBehavior}

----

<sup>

1. For this reason, linear B-splines are often referred to as *2nd order*, which may be a source of confusion since the interpolating function itself is linear, i.e. of first order. In `Interpolations.jl`, we will try to avoid this confusion by referring to interpolation degree by "linear", "quadratic" etc.

2. Although the separation of concepts 1-3 from extrapolation behavior is computationally sound, it does allow for some interesting, yet probably nonsensical, combinations of interpolation degrees, boundary conditions and extrapolation behavior. One could imagine for example a constant interpolation (which needs no boundary condition) with linear extrapolation, in which case the interpolating function is a sequence of "steps" with 0-derivative inside the domain, while suddenly having a nonzero derivative outside. Similarly, in most cases where the extrapolation behavior is defined as constant or reflecting, it will make sense to specify matching boundary conditions, but other combinations are entirely supported by `Interpolations.jl`; your milage may vary, very much...

</sup>
