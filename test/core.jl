using Interpolations, Test

@testset "==" begin

    # issue #333
    knots = ([1,1.2],)
    vals = [1.0, 2.0]
    scheme = Gridded(Linear())
    itp = interpolate(knots, vals , scheme)
    @test itp == itp
    @test itp == deepcopy(itp)
    @test itp != interpolate(([1,1.3],), vals, scheme)
    @test itp != interpolate(knots, [1.0, 3.0], scheme)
    @test itp != interpolate(knots, vals, Gridded(Constant()))

    scheme = BSpline(Quadratic(Reflect(OnCell())))
    vals = [1.0, 2.0]
    itp = interpolate(vals, scheme)
    @test itp == itp
    @test itp == deepcopy(itp)
    @test itp != interpolate(vals, BSpline(Linear()))
    @test itp != interpolate(randn(2), scheme)

    vals = rand(Float64, 2,3)
    itp1 = interpolate(vals, BSpline(Quadratic(Flat(OnGrid()))))
    itp2 = interpolate(Float64.(vals), BSpline(Quadratic(Flat(OnGrid()))))
    @test itp1 == itp2
end

@testset "Core" begin
    A = reshape([0], 1, 1, 1, 1, 1)
    IA = Interpolations.InterpGetindex(A)
    wis = ntuple(d->Interpolations.WeightedAdjIndex(1, (1,)), ndims(A))
    @test @inferred(IA[wis...]) === 0
    wis = ntuple(d->Interpolations.WeightedAdjIndex(1, (1.0,)), ndims(A))
    @test @inferred(IA[wis...]) === 0.0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,), (1,)), ndims(A))
    @test @inferred(IA[wis...]) === 0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,), (1.0,)), ndims(A))
    @test @inferred(IA[wis...]) === 0.0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,1), (1,0)), ndims(A))
    @test @inferred(IA[wis...]) === 0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,1), (1.0,0.0)), ndims(A))
    @test @inferred(IA[wis...]) === 0.0

    wi = Interpolations.WeightedAdjIndex(2, (0.2, 0.8))
    @test wi*2 === Interpolations.WeightedAdjIndex(2, (0.4, 1.6))
    @test wi/2 === Interpolations.WeightedAdjIndex(2, (0.1, 0.4))
    @test Interpolations.indextuple(wi) == (2, 3)
    wi = Interpolations.WeightedArbIndex((2, -1), (0.2, 0.8))
    @test wi*2 === Interpolations.WeightedArbIndex((2, -1), (0.4, 1.6))
    @test wi/2 === Interpolations.WeightedArbIndex((2, -1), (0.1, 0.4))
    @test Interpolations.indextuple(wi) == (2, -1)
end

@testset "Hash" begin
    # Issue 339: Base.hash accesses out-of-bound indexes due to AbstractArray interface
    # Relevant to Distributed computing
    xr = collect(0.0:1.0:5.0)
    yr = collect(1.0:1.0:6.0)
    itp = interpolate((xr,),yr,Gridded(Linear()))
    # using Distributed
    # addprocs(2)
    # @everywhere begin
    #    using Interpolations
    #    callitp(itp, x) = itp(x)
    # end
    # r = remotecall(callitp, 2, itp, 1.2)
    # fetch(r) # 2.2
    @test hash(itp) != 0
    etp = linear_interpolation([2, 3], [4, 5])
    @test hash(etp.itp) != 0
    @test hash(etp) != 0
end

@testset "BasicInferability" begin
    itp = cubic_spline_interpolation(1:3, randn(3))
    @test (@inferred Interpolations.itptype(itp)) == BSpline{Cubic{Line{OnGrid}}}
    @test (@inferred eltype(itp)) == Float64
    @test (@inferred ndims(itp)) == 1
end

@testset "checkbounds" begin
    itp = interpolate(randn(2,2,1), (BSpline(Linear()), BSpline(Linear()), NoInterp()))
    @test !Base.checkbounds(Bool, itp, 1)
    @test Base.checkbounds(Bool, itp, 1, 1)
    @test Base.checkbounds(Bool, itp, 1, 1, 1)
    @test Base.checkbounds(Bool, itp, 1, 1, 1, 1)
    @test !Base.checkbounds(Bool, itp, 1, 1, 1, 1, 2)
end
