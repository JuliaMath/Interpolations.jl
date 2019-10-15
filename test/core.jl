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
    wis = ntuple(d->Interpolations.WeightedAdjIndex(1, (1,)), ndims(A))
    @test @inferred(A[wis...]) === 0
    wis = ntuple(d->Interpolations.WeightedAdjIndex(1, (1.0,)), ndims(A))
    @test @inferred(A[wis...]) === 0.0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,), (1,)), ndims(A))
    @test @inferred(A[wis...]) === 0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,), (1.0,)), ndims(A))
    @test @inferred(A[wis...]) === 0.0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,1), (1,0)), ndims(A))
    @test @inferred(A[wis...]) === 0
    wis = ntuple(d->Interpolations.WeightedArbIndex((1,1), (1.0,0.0)), ndims(A))
    @test @inferred(A[wis...]) === 0.0

    wi = Interpolations.WeightedAdjIndex(2, (0.2, 0.8))
    @test wi*2 === Interpolations.WeightedAdjIndex(2, (0.4, 1.6))
    @test wi/2 === Interpolations.WeightedAdjIndex(2, (0.1, 0.4))
    @test Interpolations.indextuple(wi) == (2, 3)
    wi = Interpolations.WeightedArbIndex((2, -1), (0.2, 0.8))
    @test wi*2 === Interpolations.WeightedArbIndex((2, -1), (0.4, 1.6))
    @test wi/2 === Interpolations.WeightedArbIndex((2, -1), (0.1, 0.4))
    @test Interpolations.indextuple(wi) == (2, -1)
end
