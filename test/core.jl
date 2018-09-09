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
end
