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

end
