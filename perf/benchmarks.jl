using BenchmarkTools, Base.Cartesian
using Interpolations

const suite = BenchmarkGroup()

suite["bsplines"] = BenchmarkGroup()
for bspline in ["constant", "linear", "quadratic", "cubic"]
    suite["bsplines"][bspline] = BenchmarkGroup()
end

# To evaluate at fractional positions without any "unnecessary"
# overhead, safer to use Base.Cartesian
@generated function sumvalues(itp::AbstractInterpolation{T,N}, inds) where {T,N}
    quote
        @nexprs $N d->inds_d = inds[d]
        s = zero(eltype(itp))
        @inbounds @nloops $N i d->inds_d begin
            s += @ncall($N, itp, i)
        end
        s
    end
end

function sumvalues_indices(itp)
    inds = axes(itp)
    n = Int(round(10^(3/ndims(itp))))
    ntuple(d->collect(range(first(inds[d])+0.001, stop=last(inds[d])-0.001, length=n)), ndims(itp))
end

strip_prefix(str::AbstractString) = replace(str, "Interpolations."=>"")
benchstr(::Type{T}) where {T<:Interpolations.GridType} = strip_prefix(string(T))

benchstr(::Type{Constant}) = "Constant()"
benchstr(::Type{Linear}) = "Linear()"
benchstr(::Type{Quadratic{BC}}, ::Type{GT}) where {BC<:Interpolations.BoundaryCondition,GT<:Interpolations.GridType} =
    string("Quadratic(", strip_prefix(string(BC)), "(", strip_prefix(string(GT)), "()))")
benchstr(::Type{Cubic{BC}}, ::Type{GT}) where {BC<:Interpolations.BoundaryCondition,GT<:Interpolations.GridType} =
    string("Cubic(", strip_prefix(string(BC)), "(", strip_prefix(string(GT)), "()))")

groupstr(::Type{Constant}) = "constant"
groupstr(::Type{Linear}) = "linear"
groupstr(::Type{Quadratic}) = "quadratic"
groupstr(::Type{Cubic}) = "cubic"


for A in (collect(Float64, 1:3),
          reshape(collect(Float64, 1:9), 3, 3),
          reshape(collect(Float64, 1:27), 3, 3, 3))
    # Constant & Linear
    for D in (Constant, Linear)
        gstr = groupstr(D)
        Ac = copy(A)
        idstr = string(ndims(A), "d_", benchstr(D), '_', benchstr(OnGrid))
        suite["bsplines"][gstr][string(idstr, "_construct")] =
            @benchmarkable interpolate($Ac, BSpline($D()))
        itp = interpolate(copy(A), BSpline(D()))
        inds = sumvalues_indices(itp)
        suite["bsplines"][gstr][string(idstr, "_use")] =
            @benchmarkable sumvalues($itp, $inds)
    end
    # Quadratic
    gstr = groupstr(Quadratic)
    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        Ac = copy(A)
        idstr = string(ndims(A), "d_", benchstr(Quadratic{BC}, GT))
        suite["bsplines"][gstr][string(idstr, "_construct")] =
            @benchmarkable interpolate($Ac, BSpline(Quadratic($BC($GT()))))
        itp = interpolate(copy(A), BSpline(Quadratic(BC(GT()))))
        inds = sumvalues_indices(itp)
        suite["bsplines"][gstr][string(idstr, "_use")] =
            @benchmarkable sumvalues($itp, $inds)
    end
    for BC in (InPlace,InPlaceQ)
        Ac = copy(A)
        idstr = string(ndims(A), "d_", benchstr(Quadratic{BC}, OnCell))
        suite["bsplines"][gstr][string(idstr, "_construct")] =
            @benchmarkable interpolate!($Ac, BSpline(Quadratic($BC(OnCell()))))
        itp = interpolate!(copy(A), BSpline(Quadratic(BC(OnCell()))))
        inds = sumvalues_indices(itp)
        suite["bsplines"][gstr][string(idstr, "_use")] =
            @benchmarkable sumvalues($itp, $inds)
    end
    # Cubic
    gstr = groupstr(Cubic)
    for BC in (Flat,Line,Free,Periodic), GT in (OnGrid, OnCell)
        Ac = copy(A)
        idstr = string(ndims(A), "d_", benchstr(Cubic{BC}, GT))
        suite["bsplines"][gstr][string(idstr, "_construct")] =
            @benchmarkable interpolate($Ac, BSpline(Cubic($BC($GT()))))
        itp = interpolate(copy(A), BSpline(Cubic(BC(GT()))))
        inds = sumvalues_indices(itp)
        suite["bsplines"][gstr][string(idstr, "_use")] =
            @benchmarkable sumvalues($itp, $inds)
    end
end

paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals);
else
    @info "Tuning suite (this may take a while)"
    tune!(suite)
    BenchmarkTools.save(paramspath, params(suite));
end

# To run the benchmarks:
# results = run(suite, verbose = true, seconds = 1)
# BenchmarkTools.save(filename, results)
