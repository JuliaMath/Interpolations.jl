using Test
using Interpolations

@testset "iterate - 1d interpolation" begin
    itp = interpolate(1:5, BSpline(Linear()))
    kiter = knots(itp)

    @test Base.IteratorSize(kiter) === Base.HasLength()
    @test length(kiter) == 5
    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) <: Int
    @test collect(kiter) == [1, 2, 3, 4, 5]

    # Sample for loop
    kout = []
    for (k,) ∈ kiter
        push!(kout, k)
    end
    @test kout == collect(kiter)
end

@testset "iterate - 2d interpolation" begin
    itp = interpolate(rand(3,3), BSpline(Constant()))
    kiter = knots(itp)

    @test Base.IteratorSize(kiter) == Base.HasShape{2}()
    @test length(kiter) == 9
    @test size(kiter) == (3, 3)
    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) <: Tuple{Int, Int}

    k = collect(kiter)
    @test length(k) == 9
    @test size(k) == (3,3)
    @test eltype(k) <: Tuple{Int, Int}
    @test k == Iterators.product(1:3, 1:3) |> collect

    # Sample for loop
    kout = []
    for (kx, ky) ∈ kiter
        push!(kout, (kx, ky))
    end
    @test kout == Iterators.product(1:3, 1:3) |> collect |> vec

end


@testset "1D - iterate - $itp" for itp ∈ [ BSpline(Constant()) ]
    x = 1.0:5.0
    itp = interpolate(x.^2, itp)
    @test typeof(itp) <: AbstractInterpolation
    kiter = knots(itp)
    @test eltype(kiter) <: Int
    @test length(kiter) == 5
    @test collect(kiter) == x

    # Non repeating knots
    @testset "extrapolation - $etp" for etp ∈ [ Throw(), Flat(), Line()]
        extrp = extrapolate(itp, etp)
        kiter = knots(extrp)
        @test eltype(kiter) <: Int
        @test collect(kiter) == x
    end

    @testset "extrapolation - Periodic" begin
        extrp = extrapolate(itp, Periodic())
        k = knots(extrp)
        k10 = Iterators.take(k, 10) |> collect
        @test k10 ≈ collect(1.0:10.0)
        @test extrp.(k10) ≈ vcat(x[1:end-1].^2, x[1:end-1].^2, x[1:2].^2)
    end
    @testset "extrapolation - Reflect" begin
        extrp = extrapolate(itp, Reflect())
        k = knots(extrp)
        k10 = Iterators.take(k, 10) |> collect
        @test k10 ≈ collect(1.0:10.0)
        @test extrp.(k10) ≈ vcat(x.^2, reverse(x[1:end - 1].^2), x[2].^2)
    end
end

@testset "2D - iteration - Bounded" begin
    etp = LinearInterpolation(([1, 2, 3], [1, 2, 3]), rand(3, 3);
        extrapolation_bc=(Throw(), Line())
    )
    kiter = knots(etp)
    @test Base.IteratorSize(kiter) == Base.HasShape{2}()
    @test length(kiter) == 9
    @test size(kiter) == (3, 3)
    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) == Tuple{Int,Int}

    k = Iterators.take(kiter, 5) |> collect
    @test length(k) == 5
    @test typeof(k) <: Vector{Tuple{Int,Int}}
    kexp = Iterators.product(1:3, 1:3) |> x -> Iterators.take(x, 5) |> collect
    @test k == kexp
end

@testset "2D - iteration - Unbounded" begin
    etp = LinearInterpolation(([1, 2, 3], [1, 2, 3]), rand(3, 3);
        extrapolation_bc=(Line(), Periodic(OnGrid()))
    )
    kiter = knots(etp)
    @test Base.IteratorSize(kiter) == Base.IsInfinite()
    @test_throws ArgumentError size(kiter)
    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) == Tuple{Int,Int}

    k = Iterators.take(kiter, 20) |> collect
    @test length(k) == 20
    @test typeof(k) <: Vector{Tuple{Int,Int}}
    kexp = Iterators.product(1:3, Iterators.countfrom(1)) |>
        x -> Iterators.take(x, 20) |> collect
    @test k == kexp
end
