using Interpolations: degree,
                      getknots,
                      coefficients,
                      itpflag,
                      DimSpec,
                      AbstractLanczos,
                      lbounds,
                      ubounds

@testset "Lanczos" begin

@testset "Lanczos($N)" for N in 2:4
    X = 1:100
    X = [X; reverse(X)[2:end]]
    itp = interpolate(X, Lanczos(N))

    # properties
    @test getknots(itp) == axes(itp) == axes(X)
    @test coefficients(itp)[axes(itp)...] == X
    @test itpflag(itp) isa Lanczos{N}
    @test size(itp) == size(X)
    @test lbounds(itp) == (firstindex(X),)
    @test ubounds(itp) == (lastindex(X),)

    # recovers input points exactly
    @test itp.(X) == X

    @test 1 < itp(1.5) < 2
    @test 99 < itp(99.5) < 100

    # symmetry check
    interpolant = itp.(X)
    @test interpolant ≈ reverse(interpolant)
end

@testset "Lanczos OpenCV4" begin
    X = 1:100
    X = [X; reverse(X)[2:end]]
    itp = interpolate(X, Lanczos4OpenCV())

    # properties
    @test getknots(itp) == axes(itp) == axes(X)
    @test coefficients(itp)[axes(itp)...] == X
    @test itpflag(itp) isa Lanczos4OpenCV
    @test size(itp) == size(X)
    @test lbounds(itp) == (firstindex(X),)
    @test ubounds(itp) == (lastindex(X),)

    # recovers input points exactly
    @test itp.(X) == X

    @test 1 < itp(1.5) < 2
    @test 99 < itp(99.5) < 100
    # Test precision of Lanczos OpenCV4 #451, test not NaN
    @test itp(nextfloat(1.0)) ≈ itp(1.0)


    # symmetry check
    interpolant = itp.(X)
    @test interpolant ≈ reverse(interpolant)
end

@testset "Lanczos OpenCV4 Faithful" begin
    idx = Interpolations.degree(Lanczos4OpenCVFaithful{Float32}())
    @test 1f0 == Interpolations.value_weights(Lanczos4OpenCVFaithful{Float32}(), 0.0f0)[idx]
    idx = Interpolations.degree(Lanczos4OpenCVFaithful{Float64}())
    @test 1.0 == Interpolations.value_weights(Lanczos4OpenCVFaithful{Float64}(), 0.0)[idx]
    idx = Interpolations.degree(Lanczos4OpenCV())
    @test 1.0 == Interpolations.value_weights(Lanczos4OpenCV(), 0.0)[idx]
    # were _lanczos4_opencv to be refactored to use Float32 throughout,
    # the 4th element returned would be 0.0
end

end
