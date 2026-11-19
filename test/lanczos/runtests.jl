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

@testset "Scaled Lanczos bounds" begin
    signals = (ones(11), collect(0.0:10.0), [sin(i) for i in 1:11])
    physical_ranges = (0:0.1:1, range(-2.0, 4.0; length=11), 10:2:30)
    for kernel in (Lanczos(2), Lanczos(4), Lanczos4OpenCV())
        for values in signals, r in physical_ranges
            unscaled = interpolate(values, kernel)
            scaled = scale(unscaled, r)
            @test lbounds(scaled) == (first(r),)
            @test ubounds(scaled) == (last(r),)
            @test bounds(scaled) == ((first(r), last(r)),)
            @test scaled.(r) ≈ values
            for x in range(first(r), last(r); length=23)
                index = x == first(r) ? firstindex(values) :
                        x == last(r) ? lastindex(values) :
                        (x - first(r)) / step(r) + firstindex(values)
                @test scaled(x) ≈ unscaled(index)
            end
            @test_throws BoundsError scaled(prevfloat(Float64(first(r))))
            @test_throws BoundsError scaled(nextfloat(Float64(last(r))))
            @test lbounds(unscaled) == (firstindex(values),)
            @test ubounds(unscaled) == (lastindex(values),)
            @test_throws BoundsError unscaled(0.0)
            @test_throws BoundsError unscaled(12.0)
            @test_throws ArgumentError scale(unscaled, reverse(r))
            @test_throws ArgumentError scale(unscaled, range(first(r), last(r); length=10))
        end

        values = reshape(collect(1.0:20.0), 4, 5)
        xs, ys = range(-2.0, 4.0; length=4), range(10.0, 20.0; length=5)
        unscaled = interpolate(values, kernel)
        scaled = scale(unscaled, xs, ys)
        @test bounds(scaled) == ((first(xs), last(xs)), (first(ys), last(ys)))
        @test scaled(xs, ys) ≈ values
        for x in (first(xs), 0.33, last(xs)), y in (first(ys), 13.3, last(ys))
            i = (x - first(xs)) / step(xs) + firstindex(values, 1)
            j = (y - first(ys)) / step(ys) + firstindex(values, 2)
            @test scaled(x, y) ≈ unscaled(i, j)
        end
        @test_throws BoundsError scaled(prevfloat(first(xs)), first(ys))
        @test_throws BoundsError scaled(first(xs), nextfloat(last(ys)))
    end
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
