@testset "Multivalued" begin
    # 1d
    A0 = rand(20)
    A = reinterpret(MyPair{Float64}, A0)
    a1, a2 = A0[1:2:end], A0[2:2:end]
    @test length(A) == 10
    itp = interpolate(A, BSpline(Constant()))
    @test itp(3.2) ≈ MyPair(A0[5],A0[6])
    itp = interpolate(A, BSpline(Linear()))
    @test itp(3.2) ≈ 0.8*MyPair(A0[5],A0[6]) + 0.2*MyPair(A0[7],A0[8])
    it = BSpline(Quadratic(Flat(OnGrid())))
    itp = interpolate(A, it)
    @test itp(3.2) ≈ MyPair(interpolate(a1, it)(3.2), interpolate(a2, it)(3.2))

    # 2d
    A0 = rand(100)
    A = reshape(reinterpret(MyPair{Float64}, A0), (10,5))
    a1, a2 = reshape(A0[1:2:end], (10,5)), reshape(A0[2:2:end], (10,5))
    for it in (BSpline(Constant()),
               BSpline(Linear()),
               BSpline(Quadratic(Flat(OnGrid()))))
        itp = interpolate(A, it)
        @test itp(3.2,1.8) ≈ MyPair(interpolate(a1, it)(3.2,1.8), interpolate(a2, it)(3.2,1.8))
    end
end
