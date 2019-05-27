using Interpolations, Test

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
