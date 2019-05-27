import Interpolations, Grid, Dierckx, GridInterpolations, ApproXD
using PyCall
@pyimport scipy.interpolate as py
include(joinpath(JULIA_HOME, "..", "..", "examples", "ndgrid.jl"))

# For the evaluation points, loop over 2:size(itp,d)-1.
# Since some interpolants are fast, we must ensure this call is type-stable.
@generated function iterrange(itp::AbstractArray{T,N}) where {T,N}
    min = ntuple(d->2, N)
    maxex = Expr(:tuple, Expr[:(size(itp,$d)-1) for d = 1:N]...)
    :(CartesianRange(CartesianIndex($min), CartesianIndex($maxex)))
end

# These have a function-barrier so are less critical
make_knots(A) = ntuple(d->collect(linspace(1,size(A,d),size(A,d))), ndims(A))
make_xi(A)    = ntuple(d->collect(linspace(2,size(A,d)-1,size(A,d)-2)), ndims(A))

## Interpolations and Grid
function evaluate_grid(itp::Union{Array,Interpolations.AbstractInterpolation,Grid.InterpGrid}, A)
    s = zero(eltype(itp)) + zero(eltype(itp))
    for I in iterrange(itp)
        s += itp[I]
    end
    s
end

# Fast approach for GriddedInterpolation
function evaluate_grid(itp::Interpolations.GriddedInterpolation, A)
    V = itp[make_xi(A)...]
    sum(V)
end

# Slow approach for GriddedInterpolation
function evaluate_grid_scalar(itp::Interpolations.GriddedInterpolation, A)
    s = zero(eltype(itp)) + zero(eltype(itp))
    for I in iterrange(itp)
        s += itp[I]
    end
    s
end

## Dierckx
dierckx_constr(A::AbstractVector, k) = Dierckx.Spline1D(collect(1:length(A)), A, k=k)
dierckx_constr(A::Matrix, k) = Dierckx.Spline2D(collect(1:size(A,1)), collect(1:size(A,2)), A, kx=k, ky=k)

Dierckx.evaluate(itp::Dierckx.Spline1D, I::CartesianIndex{1}) = Dierckx.evaluate(itp, I[1])
Dierckx.evaluate(itp::Dierckx.Spline2D, I::CartesianIndex{2}) = Dierckx.evaluate(itp, I[1], I[2])

# Fast approach for Dierckx
function evaluate_grid(itp::Dierckx.Spline1D, A)
    V = Dierckx.evaluate(itp, make_xi(A)...)
    sum(V)
end

function evaluate_grid(itp::Dierckx.Spline2D, A)
    V = Dierckx.evalgrid(itp, make_xi(A)...)
    sum(V)
end

# Slow approach for Dierckx
function evaluate_grid_scalar(itp::Union{Dierckx.Spline1D,Dierckx.Spline2D}, A)
    T = eltype(A)
    s = zero(T) + zero(T)
    for I in iterrange(A)
        s += Dierckx.evaluate(itp, I)
    end
    s
end

## GridInterpolations
gi_constr(A) = GridInterpolations.RectangleGrid(make_knots(A)...)

stuff!(x, index::CartesianIndex{1}) = (x[1] = index[1])
stuff!(x, index::CartesianIndex{2}) = (x[1] = index[1]; x[2] = index[2])
stuff!(x, index::CartesianIndex{3}) = (x[1] = index[1]; x[2] = index[2]; x[3] = index[3])
stuff!(x, index::CartesianIndex{4}) = (x[1] = index[1]; x[2] = index[2]; x[3] = index[3]; x[4] = index[4])

function evaluate_grid(grid::GridInterpolations.RectangleGrid, A)
    T = eltype(A)
    x = Array{eltype(T)}( ndims(A))  # in case T is RGB{Float32}
    s = zero(T) + zero(T)
    vA = vec(A)
    for I in iterrange(A)
        stuff!(x, I)
        s += GridInterpolations.interpolate(grid, vA, x)
    end
    s
end

## ApproXD
ax_constr(A) = ApproXD.Lininterp(A, Vector{Float64}[make_knots(A)...])

function evaluate_grid(grid::ApproXD.Lininterp, A)
    T = eltype(A)
    x = Array{eltype(T)}( ndims(A))  # in case T is RGB{Float32}
    s = zero(T) + zero(T)
    vA = vec(A)
    result = Array{T}(1)
    which = [1]
    for I in iterrange(A)
        stuff!(x, I)
        ApproXD.getValue!(result, grid, x, which)
        s += result[1]
    end
    s
end

## Python's RegularGridInterpolator
pyindexes(A::Vector) = collect(2:length(A)-1)
function pyindexes(A)
    g = ndgrid(make_xi(A)...)
    hcat(map(vec, g)...)
end

function evaluate_grid(itp::PyObject, A)
    indexes = pyindexes(A)
    V = pycall(itp, PyAny, indexes)
    sum(V)
end

## Error-handling
function msg(err)
    if isdefined(err, :msg)
        return err.msg
    end
    string(typeof(err))
end


## OK, let's do it!
const npoints = 10^6
const maxdims = 4
# const eltypes = (Float64,Float32, RGB{Float32})
const eltypes = (Float64,)
constr_eval_name_ord = ((A -> Interpolations.interpolate(A, Interpolations.BSpline(Interpolations.Constant()), Interpolations.OnCell()),evaluate_grid,"IBSpline",0),
                    (A -> Interpolations.interpolate(A, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid()),evaluate_grid,"IBSpline",1),
                    (A -> Interpolations.interpolate(A, Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Flat())), Interpolations.OnCell()),evaluate_grid,"IBSpline",2),
                    (A -> Grid.InterpGrid(A, Grid.BCnil, Grid.InterpNearest),evaluate_grid,"Grid",0),
                    (A -> Grid.InterpGrid(A, Grid.BCnil, Grid.InterpLinear),evaluate_grid,"Grid",1),
                    (A -> Grid.InterpGrid(A, Grid.BCnil, Grid.InterpQuadratic),evaluate_grid,"Grid",2),
                    (A -> Interpolations.interpolate(make_knots(A), A, Interpolations.Gridded(Interpolations.Constant())),evaluate_grid,"IGridded vec",0),
                    (A -> Interpolations.interpolate(make_knots(A), A, Interpolations.Gridded(Interpolations.Constant())),evaluate_grid_scalar,"IGridded scalar",0),
                    (A -> Interpolations.interpolate(make_knots(A), A, Interpolations.Gridded(Interpolations.Linear())),evaluate_grid,"IGridded vec",1),
                    (A -> Interpolations.interpolate(make_knots(A), A, Interpolations.Gridded(Interpolations.Linear())),evaluate_grid_scalar,"IGridded scalar",1),
                    (A -> dierckx_constr(A, 0),evaluate_grid,"Dierckx vec",0),
                    (A -> dierckx_constr(A, 1),evaluate_grid,"Dierckx vec",1),
                    (A -> dierckx_constr(A, 1),evaluate_grid_scalar,"Dierckx scalar",1),
                    (A -> dierckx_constr(A, 2),evaluate_grid,"Dierckx vec",2),
                    (A -> dierckx_constr(A, 2),evaluate_grid_scalar,"Dierckx scalar",2),
                    (A -> gi_constr(A),evaluate_grid,"GI",1),
                    (A -> ax_constr(A),evaluate_grid,"ApproXD",1),
                    (A -> py.RegularGridInterpolator(make_knots(A), A, method="nearest"),evaluate_grid,"Py vec",0),
                    (A -> py.RegularGridInterpolator(make_knots(A), A, method="linear"),evaluate_grid,"Py vec",1),
                    )
testnames = [c_e_v[3] for c_e_v in constr_eval_name_ord]
testord   = [c_e_v[4] for c_e_v in constr_eval_name_ord]
tconstr = fill(NaN, length(constr_eval_name_ord), maxdims, length(eltypes))
teval   = fill(NaN, length(constr_eval_name_ord), maxdims, length(eltypes))
for ieltype = 1:length(eltypes)
    etyp = eltypes[ieltype]
    println(etyp)
    for ndim = 1:maxdims
        println("  dimension ", ndim)
        n = round(Int, npoints^(1/ndim))
        sz = fill(n, ndim)
        A = rand(etyp, sz...)
        starget = evaluate_grid(A, A)
        for (i,c_e_n) in enumerate(constr_eval_name_ord)
            fc, ev, name, ord = c_e_n
            println("    ", name, " ", ord)
            try
                itp = fc(A)   # JIT-compile
                gc()
                tic(); fc(A); tconstr[i, ndim, ieltype] = toq()
                s = ev(itp, A)  # JIT-compile
                @assert abs(s - starget) < 1000*eps(abs(starget))
                gc()
                tic(); ev(itp, A); teval[i, ndim, ieltype] = toq()
            catch err
                warn("skipping because of error: ", msg(err))
            end
        end
    end
end