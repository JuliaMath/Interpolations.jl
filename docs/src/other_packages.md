# Other Interpolation Packages

Other interpolation packages for Julia include:

- [ApproXD.jl](https://github.com/floswald/ApproXD.jl) implements B-spline and linear interpolation in Julia.
- [BarycentricInterpolation.jl](https://github.com/dawbarton/BarycentricInterpolation.jl) implements the Barycentric formula for polynomial interpolation on equispaced points and Chebyshev points of the first and second kind.
- [BasicInterpolators.jl](https://github.com/markmbaum/BasicInterpolators.jl) provides a collection of common interpolation recipes for basic applications.
- [BSplineKit.jl](https://github.com/jipolanco/BSplineKit.jl) offers tools for B-spline based Galerkin and collocation methods, including for interpolation and approximation.
- [Curves.jl](https://github.com/lungben/Curves.jl) supports log-interpolation via immutable `Curve` objects.
- [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) is a library for performing interpolations of one-dimensional data.
- [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl) is a wrapper for the dierckx Fortran library, which also underlies `scipy.interpolate`.
- [DIVAnd.jl](https://github.com/gher-ulg/DIVAnd.jl) for N-dimensional smoothing interpolation. 
- [FastChebInterp.jl](https://github.com/stevengj/FastChebInterp.jl) does fast multidimensional Chebyshev interpolation on a hypercube using separable grid of interpolation points.
- [FEMBasis.jl](https://github.com/JuliaFEM/FEMBasis.jl) contains interpolation routines for standard finite element function spaces.
- [FineShift.jl](https://github.com/emmt/FineShift.jl) does fast sub-sample shifting of multidimensional arrays.
- [FourierTools.jl](https://github.com/bionanoimaging/FourierTools.jl) includes sinc interpolation for up and down sampling.
- [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) provides interpolation and simulation methods over complex 2D and 3D meshes.
- [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl) performs multivariate interpolation on a rectilinear grid.
- [InterpolationKernels.jl](https://github.com/emmt/InterpolationKernels.jl) provides a library of interpolation kernels.
- [KernelInterpolation.jl](https://github.com/JoshuaLampert/KernelInterpolation.jl) implements scattered data interpolations in arbitrary dimensions by radial basis functions with support for solving linear partial differential equations.
- [KissSmoothing.jl](https://github.com/francescoalemanno/KissSmoothing.jl) implements denoising and a Radial Basis Function estimation procedure.
- [LinearInterpolations.jl](https://github.com/jw3126/LinearInterpolations.jl) allows for interpolation using weighted averages allowing probability distributions, rotations, and other Lie groups to be interpolated.
- [LinearInterpolators.jl](https://github.com/emmt/LinearInterpolators.jl) provides linear interpolation methods for Julia based on InterpolationKernels.jl, above.
- [LocalFunctionApproximation.jl](https://github.com/sisl/LocalFunctionApproximation.jl) provides local function approximators that interpolates a scalar-valued function across a vector space.
- [NaturalNeighbours.jl](https://github.com/DanielVandH/NaturalNeighbours.jl) provides natural neighbour interpolation methods for scattered two-dimensional point sets, with support for derivative generation.
- [PCHIPInterpolation.jl](https://github.com/gerlero/PCHIPInterpolation.jl) for monotonic interpolation.
- [PiecewiseLinearApprox.jl](https://github.com/RJDennis/PiecewiseLinearApprox.jl) performs piecewise linear interpolation over an arbitrary number of dimensions.
- [ScatteredInterpolation.jl](https://github.com/eljungsk/ScatteredInterpolation.jl) interpolates scattered data in arbitrary dimensions.

Some of these packages support methods that `Interpolations` does not,
so if you can't find what you need here, check one of them or submit a
pull request here.

If you would like to list a registered package that is related to interpolation, please create a Github issue.


