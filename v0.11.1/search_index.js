var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Interpolations-1",
    "page": "Home",
    "title": "Interpolations",
    "category": "section",
    "text": "(Image: Build Status) (Image: PkgEval Status) (Image: Interpolations)NEWS v0.9 was a breaking release. See the news for details on how to update.This package implements a variety of interpolation schemes for the Julia language.  It has the goals of ease-of-use, broad algorithmic support, and exceptional performance.Currently this package\'s support is best for B-splines and also supports irregular grids.  However, the API has been designed with intent to support more options. Pull-requests are more than welcome! It should be noted that the API may continue to evolve over time.Other interpolation packages for Julia include:Dierckx.jl\nGridInterpolations.jl\nApproXD.jlSome of these packages support methods that Interpolations does not, so if you can\'t find what you need here, check one of them or submit a pull request here."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Justusing Pkg\nPkg.add(\"Interpolations\")from the Julia REPL."
},

{
    "location": "interpolations/#",
    "page": "General usage",
    "title": "General usage",
    "category": "page",
    "text": ""
},

{
    "location": "interpolations/#General-usage-1",
    "page": "General usage",
    "title": "General usage",
    "category": "section",
    "text": "Note: the current version of Interpolations supports interpolation evaluation using index calls [], but this feature will be deprecated in future. We highly recommend function calls with () as follows.Given an AbstractArray A, construct an \"interpolation object\" itp asitp = interpolate(A, options...)where options... (discussed below) controls the type of interpolation you want to perform.  This syntax assumes that the samples in A are equally-spaced.To evaluate the interpolation at position (x, y, ...), simply dov = itp(x, y, ...)Some interpolation objects support computation of the gradient, which can be obtained asg = Interpolations.gradient(itp, x, y, ...)or asInterpolations.gradient!(g, itp, x, y, ...)where g is a pre-allocated vector.Some interpolation objects support computation of the hessian, which can be obtained ash = Interpolations.hessian(itp, x, y, ...)orInterpolations.hessian!(h, itp, x, y, ...)where h is a pre-allocated matrix.A may have any element type that supports the operations of addition and multiplication.  Examples include scalars like Float64, Int, and Rational, but also multi-valued types like RGB color vectors.Positions (x, y, ...) are n-tuples of numbers. Typically these will be real-valued (not necessarily integer-valued), but can also be of types such as DualNumbers if you want to verify the computed value of gradients. (Alternatively, verify gradients using ForwardDiff.) You can also use Julia\'s iterator objects, e.g.,function ongrid!(dest, itp)\n    for I in CartesianIndices(itp)\n        dest[I] = itp(I)\n    end\nendwould store the on-grid value at each grid point of itp in the output dest. Finally, courtesy of Julia\'s indexing rules, you can also usefine = itp(range(1,stop=10,length=1001), range(1,stop=15,length=201))There is also an abbreviated Convenience notation."
},

{
    "location": "control/#",
    "page": "Interpolation algorithms",
    "title": "Interpolation algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "control/#Control-of-interpolation-algorithm-1",
    "page": "Interpolation algorithms",
    "title": "Control of interpolation algorithm",
    "category": "section",
    "text": ""
},

{
    "location": "control/#BSplines-1",
    "page": "Interpolation algorithms",
    "title": "BSplines",
    "category": "section",
    "text": "The interpolation type is described in terms of degree and, if necessary, boundary conditions. There are currently four degrees available: Constant, Linear, Quadratic,  and Cubic corresponding to B-splines of degree 0, 1, 2, and 3 respectively.B-splines of quadratic or higher degree require solving an equation system to obtain the interpolation coefficients, and for that you must specify a boundary condition that is applied to close the system. The following boundary conditions are implemented: Flat, Line (alternatively, Natural), Free, Periodic and Reflect; their mathematical implications are described in detail in the pdf document under /doc/latex. When specifying these boundary conditions you also have to specify whether they apply at the edge grid point (OnGrid()) or beyond the edge point halfway to the next (fictitious) grid point (OnCell()).Some examples:# Nearest-neighbor interpolation\nitp = interpolate(a, BSpline(Constant()))\nv = itp(5.4)   # returns a[5]\n\n# (Multi)linear interpolation\nitp = interpolate(A, BSpline(Linear()))\nv = itp(3.2, 4.1)  # returns 0.9*(0.8*A[3,4]+0.2*A[4,4]) + 0.1*(0.8*A[3,5]+0.2*A[4,5])\n\n# Quadratic interpolation with reflecting boundary conditions\n# Quadratic is the lowest order that has continuous gradient\nitp = interpolate(A, BSpline(Quadratic(Reflect(OnCell()))))\n\n# Linear interpolation in the first dimension, and no interpolation (just lookup) in the second\nitp = interpolate(A, (BSpline(Linear()), NoInterp()))\nv = itp(3.65, 5)  # returns  0.35*A[3,5] + 0.65*A[4,5]There are more options available, for example:# In-place interpolation\nitp = interpolate!(A, BSpline(Quadratic(InPlace(OnCell()))))which destroys the input A but also does not need to allocate as much memory."
},

{
    "location": "control/#Scaled-BSplines-1",
    "page": "Interpolation algorithms",
    "title": "Scaled BSplines",
    "category": "section",
    "text": "BSplines assume your data is uniformly spaced on the grid 1:N, or its multidimensional equivalent. If you have data of the form [f(x) for x in A], you need to tell Interpolations about the grid A. If A is not uniformly spaced, you must use gridded interpolation described below. However, if A is a collection of ranges or linspaces, you can use scaled BSplines. This is more efficient because the gridded algorithm does not exploit the uniform spacing. Scaled BSplines can also be used with any spline degree available for BSplines, while gridded interpolation does not currently support quadratic or cubic splines.Some examples,A_x = 1.:2.:40.\nA = [log(x) for x in A_x]\nitp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))\nsitp = scale(itp, A_x)\nsitp(3.) # exactly log(3.)\nsitp(3.5) # approximately log(3.5)For multidimensional uniformly spaced gridsA_x1 = 1:.1:10\nA_x2 = 1:.5:20\nf(x1, x2) = log(x1+x2)\nA = [f(x1,x2) for x1 in A_x1, x2 in A_x2]\nitp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))\nsitp = scale(itp, A_x1, A_x2)\nsitp(5., 10.) # exactly log(5 + 10)\nsitp(5.6, 7.1) # approximately log(5.6 + 7.1)"
},

{
    "location": "control/#Gridded-interpolation-1",
    "page": "Interpolation algorithms",
    "title": "Gridded interpolation",
    "category": "section",
    "text": "These use a very similar syntax to BSplines, with the major exception being that one does not get to choose the grid representation (they are all OnGrid). As such one must specify a set of coordinate arrays defining the knots of the array.In 1DA = rand(20)\nA_x = collect(1.0:2.0:40.0)\nknots = (A_x,)\nitp = interpolate(knots, A, Gridded(Linear()))\nitp(2.0)The spacing between adjacent samples need not be constant, you can use the syntaxitp = interpolate(knots, A, options...)where knots = (xknots, yknots, ...) to specify the positions along each axis at which the array A is sampled for arbitrary (\"rectangular\") samplings.For example:A = rand(8,20)\nknots = ([x^2 for x = 1:8], [0.2y for y = 1:20])\nitp = interpolate(knots, A, Gridded(Linear()))\nitp(4,1.2)  # approximately A[2,6]One may also mix modes, by specifying a mode vector in the form of an explicit tuple:itp = interpolate(knots, A, (Gridded(Linear()),Gridded(Constant())))Presently there are only three modes for gridded:Gridded(Linear())whereby a linear interpolation is applied between knots,Gridded(Constant())whereby nearest neighbor interpolation is used on the applied axis,NoInterpwhereby the coordinate of the selected input vector MUST be located on a grid point. Requests for off grid coordinates results in the throwing of an error.missing data will naturally propagate through the interpolation, where some values will become missing. To avoid that, one can filter out the missing data points and use a gridded interpolation. For example:x = 1:6\nA = [i == 3 ? missing : i for i in x]\nxf = [xi for (xi,a) in zip(x, A) if !ismissing(a)]\nAf = [a for a in A if !ismissing(a)]\nitp = interpolate((xf, ), Af, Gridded(Linear()))"
},

{
    "location": "control/#Parametric-splines-1",
    "page": "Interpolation algorithms",
    "title": "Parametric splines",
    "category": "section",
    "text": "Given a set a knots with coordinates x(t) and y(t), a parametric spline S(t) = (x(t),y(t)) parametrized by t in [0,1] can be constructed with the following code adapted from a post by Tomas Lycken:using Interpolations\n\nt = 0:.1:1\nx = sin.(2π*t)\ny = cos.(2π*t)\nA = hcat(x,y)\n\nitp = scale(interpolate(A, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)\n\ntfine = 0:.01:1\nxs, ys = [itp(t,1) for t in tfine], [itp(t,2) for t in tfine]We can then plot the spline with:using Plots\n\nscatter(x, y, label=\"knots\")\nplot!(xs, ys, label=\"spline\")(Image: parametric spline)"
},

{
    "location": "extrapolation/#",
    "page": "Extrapolation",
    "title": "Extrapolation",
    "category": "page",
    "text": ""
},

{
    "location": "extrapolation/#Extrapolation-1",
    "page": "Extrapolation",
    "title": "Extrapolation",
    "category": "section",
    "text": "The call to extrapolate defines what happens if you try to index into the interpolation object with coordinates outside of its bounds in any dimension. The implemented boundary conditions are Throw, Flat, Linear, Periodic and Reflect, or you can pass a constant to be used as a \"fill\" value returned for any out-of-bounds evaluation. Periodic and Reflect require that there is a method of Base.mod that can handle the indices used.Examples:itp = interpolate(1:7, BSpline(Linear()))\netpf = extrapolate(itp, Flat())   # gives 1 on the left edge and 7 on the right edge\netp0 = extrapolate(itp, 0)        # gives 0 everywhere outside [1,7]"
},

{
    "location": "convenience-construction/#",
    "page": "Convenience Construcors",
    "title": "Convenience Construcors",
    "category": "page",
    "text": ""
},

{
    "location": "convenience-construction/#Convenience-notation-1",
    "page": "Convenience Construcors",
    "title": "Convenience notation",
    "category": "section",
    "text": "For linear and cubic spline interpolations, LinearInterpolation and CubicSplineInterpolation can be used to create interpolating and extrapolating objects handily:f(x) = log(x)\nxs = 1:0.2:5\nA = [f(x) for x in xs]\n\n# linear interpolation\ninterp_linear = LinearInterpolation(xs, A)\ninterp_linear(3) # exactly log(3)\ninterp_linear(3.1) # approximately log(3.1)\n\n# cubic spline interpolation\ninterp_cubic = CubicSplineInterpolation(xs, A)\ninterp_cubic(3) # exactly log(3)\ninterp_cubic(3.1) # approximately log(3.1)which support multidimensional data as well:f(x,y) = log(x+y)\nxs = 1:0.2:5\nys = 2:0.1:5\nA = [f(x,y) for x in xs, y in ys]\n\n# linear interpolation\ninterp_linear = LinearInterpolation((xs, ys), A)\ninterp_linear(3, 2) # exactly log(3 + 2)\ninterp_linear(3.1, 2.1) # approximately log(3.1 + 2.1)\n\n# cubic spline interpolation\ninterp_cubic = CubicSplineInterpolation((xs, ys), A)\ninterp_cubic(3, 2) # exactly log(3 + 2)\ninterp_cubic(3.1, 2.1) # approximately log(3.1 + 2.1)For extrapolation, i.e., when interpolation objects are evaluated in coordinates outside the range provided in constructors, the default option for a boundary condition is Throw so that they will return an error. Interested users can specify boundary conditions by providing an extra parameter for extrapolation_bc:f(x) = log(x)\nxs = 1:0.2:5\nA = [f(x) for x in xs]\n\n# extrapolation with linear boundary conditions\nextrap = LinearInterpolation(xs, A, extrapolation_bc = Line())\n\n@test extrap(1 - 0.2) # ≈ f(1) - (f(1.2) - f(1))\n@test extrap(5 + 0.2) # ≈ f(5) + (f(5) - f(4.8))You can also use a \"fill\" value, which gets returned whenever you ask for out-of-range values:extrap = LinearInterpolation(xs, A, extrapolation_bc = NaN)\n@test isnan(extrap(5.2))Irregular grids are supported as well; note that presently only LinearInterpolation supports irregular grids.xs = [x^2 for x = 1:0.2:5]\nA = [f(x) for x in xs]\n\n# linear interpolation\ninterp_linear = LinearInterpolation(xs, A)\ninterp_linear(1) # exactly log(1)\ninterp_linear(1.05) # approximately log(1.05)"
},

{
    "location": "api/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "api/#Interpolations.CubicSplineInterpolation",
    "page": "Library",
    "title": "Interpolations.CubicSplineInterpolation",
    "category": "function",
    "text": "etp = CubicSplineInterpolation(knots, A; bc=Line(OnGrid()), extrapolation_bc=Throw())\n\nA shorthand for extrapolate(interpolate(knots, A, BSpline(Cubic(bc))), extrapolation_bc).\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.LinearInterpolation",
    "page": "Library",
    "title": "Interpolations.LinearInterpolation",
    "category": "function",
    "text": "etp = LinearInterpolation(knots, A; extrapolation_bc=Throw())\n\nA shorthand for extrapolate(interpolate(knots, A, scheme), extrapolation_bc), where scheme is either BSpline(Linear()) or Gridded(Linear()) depending on whether knots are ranges or vectors.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.extrapolate-Union{Tuple{ET}, Tuple{IT}, Tuple{N}, Tuple{T}, Tuple{AbstractInterpolation{T,N,IT},ET}} where ET<:Union{Tuple{Vararg{Union{Tuple{BoundaryCondition,BoundaryCondition}, BoundaryCondition},N} where N}, BoundaryCondition} where IT where N where T",
    "page": "Library",
    "title": "Interpolations.extrapolate",
    "category": "method",
    "text": "extrapolate(itp, scheme) adds extrapolation behavior to an interpolation object, according to the provided scheme.\n\nThe scheme can take any of these values:\n\nThrow - throws a BoundsError for out-of-bounds indices\nFlat - for constant extrapolation, taking the closest in-bounds value\n`Line - linear extrapolation (the wrapped interpolation object must support gradient)\nReflect - reflecting extrapolation (indices must support mod)\nPeriodic - periodic extrapolation (indices must support mod)\n\nYou can also combine schemes in tuples. For example, the scheme (Line), Flat()) will use linear extrapolation in the first dimension, and constant in the second.\n\nFinally, you can specify different extrapolation behavior in different direction. ((Line),Flat()), Flat()) will extrapolate linearly in the first dimension if the index is too small, but use constant etrapolation if it is too large, and always use constant extrapolation in the second dimension.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.extrapolate-Union{Tuple{IT}, Tuple{N}, Tuple{T}, Tuple{AbstractInterpolation{T,N,IT},Any}} where IT where N where T",
    "page": "Library",
    "title": "Interpolations.extrapolate",
    "category": "method",
    "text": "extrapolate(itp, fillvalue) creates an extrapolation object that returns the fillvalue any time the indexes in itp(x1,x2,...) are out-of-bounds.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.interpolate-Union{Tuple{IT}, Tuple{AbstractArray,IT}} where IT<:Union{NoInterp, Tuple{Vararg{Union{NoInterp, BSpline},N} where N}, BSpline}",
    "page": "Library",
    "title": "Interpolations.interpolate",
    "category": "method",
    "text": "itp = interpolate(A, interpmode, gridstyle)\n\nInterpolate an array A in the mode determined by interpmode and gridstyle. interpmode may be one of\n\nBSpline(NoInterp())\nBSpline(Linear())\nBSpline(Quadratic(BC())) (see BoundaryCondition)\nBSpline(Cubic(BC()))\n\nIt may also be a tuple of such values, if you want to use different interpolation schemes along each axis.\n\ngridstyle should be one of OnGrid() or OnCell().\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.scale-Union{Tuple{IT}, Tuple{N}, Tuple{T}, Tuple{AbstractInterpolation{T,N,IT},Vararg{AbstractRange,N}}} where IT where N where T",
    "page": "Library",
    "title": "Interpolations.scale",
    "category": "method",
    "text": "scale(itp, xs, ys, ...) scales an existing interpolation object to allow for indexing using other coordinate axes than unit ranges, by wrapping the interpolation object and transforming the indices from the provided axes onto unit ranges upon indexing.\n\nThe parameters xs etc must be either ranges or linspaces, and there must be one coordinate range/linspace for each dimension of the interpolation object.\n\nFor every NoInterp dimension of the interpolation object, the range must be exactly 1:size(itp, d).\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.CardinalMonotonicInterpolation",
    "page": "Library",
    "title": "Interpolations.CardinalMonotonicInterpolation",
    "category": "type",
    "text": "CardinalMonotonicInterpolation(tension)\n\nCubic cardinal splines, uses tension parameter which must be between [0,1] Cubin cardinal splines can overshoot for non-monotonic data (increasing tension reduces overshoot).\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.Constant",
    "page": "Library",
    "title": "Interpolations.Constant",
    "category": "type",
    "text": "Constant b-splines are nearest-neighbor interpolations, and effectively return A[round(Int,x)] when interpolating.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.Cubic",
    "page": "Library",
    "title": "Interpolations.Cubic",
    "category": "type",
    "text": "Assuming uniform knots with spacing 1, the ith piece of cubic spline implemented here is defined as follows.\n\ny_i(x) = cm p(x-i) + c q(x-i) + cp q(1- (x-i)) + cpp p(1 - (x-i))\n\nwhere\n\np(δx) = 1/6 * (1-δx)^3\nq(δx) = 2/3 - δx^2 + 1/2 δx^3\n\nand the values cX for X ∈ {m, _, p, pp} are the pre-filtered coefficients.\n\nFor future reference, this expands out to the following polynomial:\n\ny_i(x) = 1/6 cm (1+i-x)^3 + c (2/3 - (x-i)^2 + 1/2 (x-i)^3) +\n         cp (2/3 - (1+i-x)^2 + 1/2 (1+i-x)^3) + 1/6 cpp (x-i)^3\n\nWhen we derive boundary conditions we will use derivatives y_0\'(x) and y_0\'\'(x)\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.FiniteDifferenceMonotonicInterpolation",
    "page": "Library",
    "title": "Interpolations.FiniteDifferenceMonotonicInterpolation",
    "category": "type",
    "text": "FiniteDifferenceMonotonicInterpolation\n\nClassic cubic interpolation, no tension parameter. Finite difference can overshoot for non-monotonic data.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.FritschButlandMonotonicInterpolation",
    "page": "Library",
    "title": "Interpolations.FritschButlandMonotonicInterpolation",
    "category": "type",
    "text": "FritschButlandMonotonicInterpolation\n\nMonotonic interpolation based on  Fritsch & Butland (1984), \"A Method for Constructing Local Monotone Piecewise Cubic Interpolants\", doi:10.1137/0905021.\n\nFaster than FritschCarlsonMonotonicInterpolation (only requires one pass) but somewhat higher apparent \"tension\".\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.FritschCarlsonMonotonicInterpolation",
    "page": "Library",
    "title": "Interpolations.FritschCarlsonMonotonicInterpolation",
    "category": "type",
    "text": "FritschCarlsonMonotonicInterpolation\n\nMonotonic interpolation based on Fritsch & Carlson (1980), \"Monotone Piecewise Cubic Interpolation\", doi:10.1137/0717021.\n\nTangents are first initialized, then adjusted if they are not monotonic can overshoot for non-monotonic data\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.Linear",
    "page": "Library",
    "title": "Interpolations.Linear",
    "category": "type",
    "text": "Assuming uniform knots with spacing 1, the ith piece of linear b-spline implemented here is defined as follows.\n\ny_i(x) = c p(x) + cp p(1-x)\n\nwhere\n\np(δx) = x\n\nand the values cX for X ∈ {_, p} are the coefficients.\n\nLinear b-splines are naturally interpolating, and require no prefiltering; there is therefore no need for boundary conditions to be provided.\n\nAlso, although the implementation is slightly different in order to re-use the framework built for general b-splines, the resulting interpolant is just a piecewise linear function connecting each pair of neighboring data points.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.LinearMonotonicInterpolation",
    "page": "Library",
    "title": "Interpolations.LinearMonotonicInterpolation",
    "category": "type",
    "text": "LinearMonotonicInterpolation\n\nSimple linear interpolation.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.Quadratic",
    "page": "Library",
    "title": "Interpolations.Quadratic",
    "category": "type",
    "text": "Assuming uniform knots with spacing 1, the ith piece of quadratic spline implemented here is defined as follows:\n\ny_i(x) = cm p(x-i) + c q(x) + cp p(1-(x-i))\n\nwhere\n\np(δx) = (δx - 1)^2 / 2\nq(δx) = 3/4 - δx^2\n\nand the values for cX for X ∈ {m,_,p} are the pre-filtered coefficients.\n\nFor future reference, this expands to the following polynomial:\n\ny_i(x) = cm * 1/2 * (x-i-1)^2 + c * (3/4 - x + i)^2 + cp * 1/2 * (x-i)^2\n\nWhen we derive boundary conditions we will use derivatives y_1\'(x-1) and y_1\'\'(x-1)\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.SteffenMonotonicInterpolation",
    "page": "Library",
    "title": "Interpolations.SteffenMonotonicInterpolation",
    "category": "type",
    "text": "SteffenMonotonicInterpolation\n\nMonotonic interpolation based on Steffen (1990), \"A Simple Method for Monotonic Interpolation in One Dimension\", http://adsabs.harvard.edu/abs/1990A%26A...239..443S\n\nOnly one pass, results usually between FritschCarlson and FritschButland.\n\n\n\n\n\n"
},

{
    "location": "api/#Public-API-1",
    "page": "Library",
    "title": "Public API",
    "category": "section",
    "text": "DocTestSetup= quote\nusing Interpolations\nendModules = [Interpolations]\nPrivate = false\nOrder = [:function, :type]"
},

{
    "location": "api/#Interpolations.boundstep",
    "page": "Library",
    "title": "Interpolations.boundstep",
    "category": "function",
    "text": "Returns half the width of one step of the range.\n\nThis function is used to calculate the upper and lower bounds of OnCell interpolation objects.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.count_interp_dims-Union{Tuple{Type{ITP}}, Tuple{ITP}} where ITP<:AbstractInterpolation",
    "page": "Library",
    "title": "Interpolations.count_interp_dims",
    "category": "method",
    "text": "n = count_interp_dims(ITP)\n\nCount the number of dimensions along which type ITP is interpolating. NoInterp dimensions do not contribute to the sum.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.gradient_weights",
    "page": "Library",
    "title": "Interpolations.gradient_weights",
    "category": "function",
    "text": "w = gradient_weights(degree, δx)\n\nCompute the weights for interpolation of the gradient at an offset δx from the \"base\" position. degree describes the interpolation scheme.\n\nExample\n\njulia> Interpolations.gradient_weights(Linear(), 0.2)\n(-1.0, 1.0)\n\nThis defines the gradient of a linear interpolation at 3.2 as y[4] - y[3].\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.hessian_weights",
    "page": "Library",
    "title": "Interpolations.hessian_weights",
    "category": "function",
    "text": "w = hessian_weights(degree, δx)\n\nCompute the weights for interpolation of the hessian at an offset δx from the \"base\" position. degree describes the interpolation scheme.\n\nExample\n\njulia> Interpolations.hessian_weights(Linear(), 0.2)\n(0.0, 0.0)\n\nLinear interpolation uses straight line segments, so the second derivative is zero.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.inner_system_diags",
    "page": "Library",
    "title": "Interpolations.inner_system_diags",
    "category": "function",
    "text": "dl, d, du = inner_system_diags{T,IT}(::Type{T}, n::Int, ::Type{IT})\n\nHelper function to generate the prefiltering equation system: generates the diagonals for a n-by-n tridiagonal matrix with eltype T corresponding to the interpolation type IT.\n\ndl, d, and du are intended to be used e.g. as in M = Tridiagonal(dl, d, du)\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.inner_system_diags-Union{Tuple{T}, Tuple{Type{T},Int64,Cubic}} where T",
    "page": "Library",
    "title": "Interpolations.inner_system_diags",
    "category": "method",
    "text": "Cubic: continuity in function value, first and second derivatives yields\n\n2/3 1/6\n1/6 2/3 1/6\n    1/6 2/3 1/6\n       ⋱  ⋱   ⋱\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "function",
    "text": "M, b = prefiltering_system{T,TC,GT<:GridType,D<:Degree}m(::T, ::Type{TC}, n::Int, ::Type{D}, ::Type{GT})\n\nGiven element types (T, TC) and interpolation scheme (GT, D) as well the number of rows in the data input (n), compute the system used to prefilter spline coefficients. Boundary conditions determine the values on the first and last rows.\n\nSome of these boundary conditions require that these rows have off-tridiagonal elements (e.g the [1,3] element of the matrix). To maintain the efficiency of solving tridiagonal systems, the Woodbury matrix identity is used to add additional elements off the main 3 diagonals.\n\nThe filtered coefficients are given by solving the equation system\n\nM * c = v + b\n\nwhere c are the sought coefficients, and v are the data points.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{BC}, Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Quadratic{BC}}} where BC<:Union{Flat{OnCell}, Reflect{OnCell}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Quadratic{Flat} OnCell and Quadratic{Reflect} OnCell amounts to setting y_1\'(x) = 0 at x=1/2. Applying this condition yields\n\n-cm + c = 0\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{BC}, Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Quadratic{BC}}} where BC<:Union{Flat{OnGrid}, Reflect{OnGrid}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Quadratic{Flat} OnGrid and Quadratic{Reflect} OnGrid amount to setting y_1\'(x) = 0 at x=1. Applying this condition yields\n\n-cm + cp = 0\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Cubic{#s12} where #s12<:Free}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Cubic{Free} OnGrid and Cubic{Free} OnCell amount to requiring an extra continuous derivative at the second-to-last cell boundary; this means y_1\'\'\'(2) = y_2\'\'\'(2), yielding\n\n1 cm -3 c + 3 cp -1 cpp = 0\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Cubic{#s12} where #s12<:Periodic}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Cubic{Periodic} OnGrid closes the system by looking at the coefficients themselves as periodic, yielding\n\nc0 = c(N+1)\n\nwhere N is the number of data points.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Cubic{Flat{OnCell}}}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Cubic{Flat}, OnCell amounts to setting y_1\'(x) = 0 at x = 1/2. Applying this condition yields\n\n-9/8 cm + 11/8 c - 3/8 cp + 1/8 cpp = 0\n\nor, equivalently,\n\n-9 cm + 11 c -3 cp + 1 cpp = 0\n\n(Note that we use y_1\'(x) although it is strictly not valid in this domain; if we were to use y_0\'(x) we would have to introduce new coefficients, so that would not close the system. Instead, we extend the outermost polynomial for an extra half-cell.)\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Cubic{Flat{OnGrid}}}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Cubic{Flat} OnGrid amounts to setting y_1\'(x) = 0 at x = 1. Applying this condition yields\n\n-cm + cp = 0\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Cubic{Line{OnCell}}}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Cubic{Line} OnCell amounts to setting y_1\'\'(x) = 0 at x = 1/2. Applying this condition yields\n\n3 cm -7 c + 5 cp -1 cpp = 0\n\n(Note that we use y_1\'(x) although it is strictly not valid in this domain; if we were to use y_0\'(x) we would have to introduce new coefficients, so that would not close the system. Instead, we extend the outermost polynomial for an extra half-cell.)\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Cubic{Line{OnGrid}}}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Cubic{Line} OnGrid amounts to setting y_1\'\'(x) = 0 at x = 1. Applying this condition gives:\n\n1 cm -2 c + 1 cp = 0\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Quadratic{#s12} where #s12<:Free}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Quadratic{Free} OnGrid and Quadratic{Free} OnCell amount to requiring an extra continuous derivative at the second-to-last cell boundary; this means that y_1\'\'(3/2) = y_2\'\'(3/2), yielding\n\n1 cm -3 c + 3 cp - cpp = 0\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Quadratic{#s12} where #s12<:Line}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Quadratic{Line} OnGrid and Quadratic{Line} OnCell amount to setting y_1\'\'(x) = 0 at x=1 and x=1/2 respectively. Since y_i\'\'(x) is independent of x for a quadratic b-spline, these both yield\n\n1 cm -2 c + 1 cp = 0\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.prefiltering_system-Union{Tuple{TC}, Tuple{T}, Tuple{Type{T},Type{TC},Int64,Quadratic{#s12} where #s12<:Periodic}} where TC where T",
    "page": "Library",
    "title": "Interpolations.prefiltering_system",
    "category": "method",
    "text": "Quadratic{Periodic} OnGrid and Quadratic{Periodic} OnCell close the system by looking at the coefficients themselves as periodic, yielding\n\nc0 = c(N+1)\n\nwhere N is the number of data points.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.rescale_gradient",
    "page": "Library",
    "title": "Interpolations.rescale_gradient",
    "category": "function",
    "text": "rescale_gradient(r::AbstractRange)\n\nImplements the chain rule dy/dx = dy/du * du/dx for use when calculating gradients with scaled interpolation objects.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.show_ranged-Tuple{IO,Any,Any}",
    "page": "Library",
    "title": "Interpolations.show_ranged",
    "category": "method",
    "text": "show_ranged(io, X, knots)\n\nA replacement for the default array-show for types that may not have the canonical evaluation points. rngs is the tuple of knots along each axis.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.value_weights",
    "page": "Library",
    "title": "Interpolations.value_weights",
    "category": "function",
    "text": "w = value_weights(degree, δx)\n\nCompute the weights for interpolation of the value at an offset δx from the \"base\" position. degree describes the interpolation scheme.\n\nExample\n\njulia> Interpolations.value_weights(Linear(), 0.2)\n(0.8, 0.2)\n\nThis corresponds to the fact that linear interpolation at x + 0.2 is 0.8*y[x] + 0.2*y[x+1].\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.BoundaryCondition",
    "page": "Library",
    "title": "Interpolations.BoundaryCondition",
    "category": "type",
    "text": "BoundaryCondition\n\nAn abstract type with one of the following values (see the help for each for details):\n\nThrow(gt)\nFlat(gt)\nLine(gt)\nFree(gt)\nPeriodic(gt)\nReflect(gt)\nInPlace(gt)\nInPlaceQ(gt)\n\nwhere gt is the grid type, e.g., OnGrid() or OnCell(). OnGrid means that the boundary condition \"activates\" at the first and/or last integer location within the interpolation region, OnCell means the interpolation extends a half-integer beyond the edge before activating the boundary condition.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.BoundsCheckStyle",
    "page": "Library",
    "title": "Interpolations.BoundsCheckStyle",
    "category": "type",
    "text": "BoundsCheckStyle(itp)\n\nA trait to determine dispatch of bounds-checking for itp. Can return NeedsCheck(), in which case bounds-checking is performed, or CheckWillPass() in which case the check will return true.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.MonotonicInterpolation",
    "page": "Library",
    "title": "Interpolations.MonotonicInterpolation",
    "category": "type",
    "text": "MonotonicInterpolation\n\nMonotonic interpolation up to third order represented by type, knots and coefficients.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.MonotonicInterpolationType",
    "page": "Library",
    "title": "Interpolations.MonotonicInterpolationType",
    "category": "type",
    "text": "MonotonicInterpolationType\n\nAbstract class for all types of monotonic interpolation.\n\n\n\n\n\n"
},

{
    "location": "api/#Interpolations.WeightedIndex",
    "page": "Library",
    "title": "Interpolations.WeightedIndex",
    "category": "type",
    "text": "wi = WeightedIndex(indexes, weights)\n\nConstruct a weighted index wi, which can be thought of as a generalization of an ordinary array index to the context of interpolation. For an ordinary vector a, a[i] extracts the element at index i. When interpolating, one is typically interested in a range of indexes and the output is some weighted combination of array values at these indexes. For example, for linear interpolation between i and i+1 we have\n\nret = (1-f)*a[i] + f*a[i]\n\nThis can be represented a[wi], where\n\nwi = WeightedIndex(i:i+1, (1-f, f))\n\ni.e.,\n\nret = sum(a[indexes] .* weights)\n\nLinear interpolation thus constructs weighted indices using a 2-tuple for weights and a length-2 indexes range. Higher-order interpolation would involve more positions and weights (e.g., 3-tuples for quadratic interpolation, 4-tuples for cubic).\n\nIn multiple dimensions, separable interpolation schemes are implemented in terms of multiple weighted indices, accessing A[wi1, wi2, ...] where each wi is the WeightedIndex along the corresponding dimension.\n\nFor value interpolation, weights will typically sum to 1. However, for gradient and Hessian computation this will not necessarily be true. For example, the gradient of one-dimensional linear interpolation can be represented as\n\ngwi = WeightedIndex(i:i+1, (-1, 1))\ng1 = a[gwi]\n\nFor a three-dimensional array A, one might compute ∂A/∂x₂ (the second component of the gradient) as A[wi1, gwi2, wi3], where wi1 and wi3 are \"value\" weights and gwi2 \"gradient\" weights.\n\nindexes may be supplied as a range or as a tuple of the same length as weights. The latter is applicable, e.g., for periodic boundary conditions.\n\n\n\n\n\n"
},

{
    "location": "api/#Internal-API-1",
    "page": "Library",
    "title": "Internal API",
    "category": "section",
    "text": "Modules = [Interpolations]\nPublic = false\nOrder = [:function, :type]"
},

]}
