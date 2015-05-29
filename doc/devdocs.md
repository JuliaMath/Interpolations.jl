# Developer documentation

Interpolations uses metaprogramming to improve performance; however, for people new to metaprogramming
this can can be a barrier.
Fortunately, with a few tips a lot of the mystique goes away.

## Looking under the hood

First let's create an interpolation object:

    julia> using Interpolations

    julia> A = rand(5)
    5-element Array{Float64,1}:
    0.74838
    0.995383
    0.978916
    0.134746
    0.430876

    julia> yitp = interpolate(A, BSpline(Linear), OnGrid)
    5-element Interpolations.BSplineInterpolation{Float64,1,Float64,Interpolations.BSpline{Interpolations.Linear},Interpolations.OnGrid}:
    0.74838
    0.995383
    0.978916
    0.134746
    0.430876

We can use this object to learn a lot about how Interpolations works.
For example, the key functionality provided by `yitp` is `getindex`, i.e., `itp[3.2]`.
Where is this implemented?

    julia> @which yitp[3.2]
    getindex{T,N}(itp::Interpolations.BSplineInterpolation{T,N,TCoefs,IT<:Interpolations.BSpline{D<:Interpolations.Degree{N}},GT<:Interpolations.GridType},xs::Real) at /home/tim/.julia/v0.4/Interpolations/src/b-splines/indexing.jl:42

Your specific output (and especially the line number) may differ, but the point is that you've now found out where this is implemented.
If you take a look at that function definition, you might see something like this:

    @generated function getindex{T,N}(itp::BSplineInterpolation{T,N}, xs::Real)
        if N > 1
            error("Linear indexing is not supported for interpolation objects")
        end
        getindex_impl(itp)
    end

This is a [generated function](http://docs.julialang.org/en/latest/manual/metaprogramming/#generated-functions), and you'll need to familiarize yourself with how these work.
The "interesting" part of the function is the call to `getindex_impl`; we can see the code that gets generated like this:

    julia> Interpolations.getindex_impl(typeof(yitp))
    quote  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/indexing.jl, line 7:
        @nexprs 1 (d->begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/indexing.jl, line 7:
                    x_d = xs[d]
                end) # line 11:
        begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 5:
            @nexprs 1 (d->begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 5:
                        begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 6:
                            ix_d = clamp(floor(Int,real(x_d)),1,size(itp,d) - 1) # line 7:
                            ixp_d = ix_d + 1 # line 8:
                            fx_d = x_d - ix_d
                        end
                    end)
        end # line 14:
        @nexprs 1 (d->begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 14:
                    begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 20:
                        c_d = 1 - fx_d # line 21:
                        cp_d = fx_d
                    end
                end) # line 17:
        @inbounds ret = c_1 * itp.coefs[ix_1] + cp_1 * itp.coefs[ixp_1] # line 18:
        ret
    end

You can see that this code makes use of [Base.Cartesian](http://docs.julialang.org/en/latest/devdocs/cartesian/), which you may also need to study.
However, the impact of these macros can be gleaned through `macroexpand`:

    julia> using Base.Cartesian

    julia> macroexpand(Interpolations.getindex_impl(typeof(yitp)))
    quote  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/indexing.jl, line 7:
        begin
            x_1 = xs[1]
        end # line 11:
        begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 5:
            begin
                begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 6:
                    ix_1 = clamp(floor(Int,real(x_1)),1,size(itp,1) - 1) # line 7:
                    ixp_1 = ix_1 + 1 # line 8:
                    fx_1 = x_1 - ix_1
                end
            end
        end # line 14:
        begin
            begin  # /home/tim/.julia/v0.4/Interpolations/src/b-splines/linear.jl, line 20:
                c_1 = 1 - fx_1 # line 21:
                cp_1 = fx_1
            end
        end # line 17:
        begin
            $(Expr(:boundscheck, false))
            begin
                ret = c_1 * itp.coefs[ix_1] + cp_1 * itp.coefs[ixp_1]
                $(Expr(:boundscheck, :(Base.pop)))
            end
        end # line 18:
        ret
    end

This is probably starting to look like something you can read. Briefly, what's happening is:

- `floor(Int,x_1)` gets clamped to the range `[1,size(itp,1)-1]`; this is the lower-bound integer grid point
- `ixp_1` is defined as `ix_1+1`; this is the upper-bound integer grid point. In Interpolations, `m` and `p` often mean "minus" and "plus", meaning the lower or upper grid point.
- The fractional part is stored in `fx_1`
- Position-coefficients `c_1` and `cp_1` associated with the lower and upper grid point are computed from `fx_1`
- The interpolation is performed using the position-coefficients, grid points, and data-coefficients and stored in `ret`, which is returned.

As useful exercises:

- Try creating a 2-dimensional linear interpolation object and examine the created code
- Create a `Quadratic` interpolation object and do the same

Once you've gotten this far, you probably understand quite a lot about how Interpolations works.
At this point, your best bet is to start looking into the helper functions used by `getindex_impl`;
once you learn how to define these, you should be able to fairly easily extend Interpolations to new types.
