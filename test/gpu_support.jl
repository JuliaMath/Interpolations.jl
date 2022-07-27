using JLArrays, Adapt
JLArrays.allowscalar(false)

@testset "1d GPU Interpolation" begin
    A_x = 1.:2.:40.
    A = [log(x) for x in A_x]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    jlitp = jl(itp)
    idx = range(2, 19., 101)
    jlidx = jl(collect(idx))
    @test itp.(idx) == collect(jlitp.(idx)) == collect(jlitp.(jlidx))
    @test gradient.(Ref(itp), idx) ==
          collect(gradient.(Ref(jlitp), idx)) ==
          collect(gradient.(Ref(jlitp), jlidx))

    sitp = scale(itp, A_x)
    jlsitp = jl(sitp)
    idx = range(1., 39., 99)
    jlidx = jl(collect(idx))
    @test sitp.(idx) == collect(jlsitp.(idx)) == collect(jlsitp.(jlidx))
    @test gradient.(Ref(sitp), idx) ==
          collect(gradient.(Ref(jlsitp), idx)) ==
          collect(gradient.(Ref(jlsitp), jlidx))


    esitp = extrapolate(sitp, Flat())
    jlesitp = jl(esitp)
    idx = range(-1., 41., 51)
    jlidx = jl(collect(idx))
    @test esitp.(idx) == collect(jlesitp.(idx)) == collect(jlesitp.(jlidx))
    @test gradient.(Ref(esitp), idx) ==
          collect(gradient.(Ref(jlesitp), idx)) ==
          collect(gradient.(Ref(jlesitp), jlidx))
end

@testset "2d GPU Interpolation" begin
    A_x = 1.:2.:40.
    A = [log(x + y) for x in A_x, y in 1.:2.:40.]
    itp = interpolate(A, (BSpline(Cubic(Line(OnGrid()))), BSpline(Linear())))
    jlitp = jl(itp)
    idx = range(2, 19., 101)
    jlidx = jl(collect(idx))
    @test itp.(idx, idx') == collect(jlitp.(idx, idx')) == collect(jlitp.(jlidx, jlidx'))
    @test gradient.(Ref(itp), idx, idx') ==
          collect(gradient.(Ref(jlitp), idx, idx')) ==
          collect(gradient.(Ref(jlitp), jlidx, jlidx'))
    @test hessian.(Ref(itp), idx, idx') ==
          collect(hessian.(Ref(jlitp), idx, idx')) ==
          collect(hessian.(Ref(jlitp), jlidx, jlidx'))

    sitp = scale(itp, A_x, A_x)
    jlsitp = jl(sitp)
    idx = range(1., 39., 99)
    jlidx = jl(collect(idx))
    @test sitp.(idx, idx') == collect(jlsitp.(idx, idx')) == collect(jlsitp.(jlidx, jlidx'))
    @test gradient.(Ref(sitp), idx, idx') ==
          collect(gradient.(Ref(jlsitp), idx, idx')) ==
          collect(gradient.(Ref(jlsitp), jlidx, jlidx'))
    @test hessian.(Ref(sitp), idx, idx') ==
          collect(hessian.(Ref(jlsitp), idx, idx')) ==
          collect(hessian.(Ref(jlsitp), jlidx, jlidx'))

    esitp = extrapolate(sitp, Flat())
    jlesitp = jl(esitp)
    idx = range(-1., 41., 51)
    jlidx = jl(collect(idx))
    @test esitp.(idx, idx') == collect(jlesitp.(idx, idx')) == collect(jlesitp.(jlidx, jlidx'))
    # gradient for `extrapolation` is currently broken under CUDA
    @test gradient.(Ref(esitp), idx, idx') ==
          collect(gradient.(Ref(jlesitp), idx, idx')) ==
          collect(gradient.(Ref(jlesitp), jlidx, jlidx'))
end

@testset "Lanczos on gpu" begin
    X = 1:100
    X = [X; reverse(X)[2:end]]
    for N = 2:4
        itp = interpolate(X, Lanczos(N))
        @test itp.(X) == collect(jl(itp).(jl(X)))
    end
    itp = interpolate(X, Lanczos4OpenCV())
    @test itp.(X) == collect(jl(itp).(jl(X)))
end

@testset "Gridded on gpu" begin
    itp1 = interpolate(Vector.((range(-1,1,101), range(-1,1,101))), randn(101,101), Gridded(Linear()))
    itp2 = interpolate((range(-1,1,101), range(-1,1,101)), randn(101,101), Gridded(Linear()))
    idx = range(-1, 1, 301)
    jlidx = jl(collect(idx))
    @test itp1.(idx, idx') == collect(jl(itp1).(idx, idx')) == collect(jl(itp1).(jlidx, jlidx'))
    @test itp2.(idx, idx') == collect(jl(itp2).(idx, idx')) == collect(jl(itp2).(jlidx, jlidx'))
end

@testset "eltype after adaption" begin
    A_x = 1.:2.:40.
    A = [log(x) for x in A_x]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    @test eltype(adapt(Array{Float32}, itp)) === Float32
    @test eltype(adapt(Array{Float32}, scale(itp, A_x))) === Float32
    @test eltype(adapt(Array{Float32}, extrapolate(scale(itp, A_x), Flat()))) === Float32
    itp = interpolate((range(-1,1,10), range(-1,1,10)), randn(10,10), Gridded(Linear()))
    @test eltype(adapt(Array{Float32}, itp)) === Float32
end
