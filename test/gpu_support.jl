using JLArrays, Adapt
JLArrays.allowscalar(false)

@testset "1d GPU Interpolation" begin
    A_x = 1.0:2.0:40.0
    A = [log(x) for x in A_x]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    jlitp = jl(itp)
    idx = 2.0:0.17:19.0
    jlidx = jl(collect(idx))
    @test itp.(idx) == collect(jlitp.(idx)) == collect(jlitp.(jlidx))
    @test Interpolations.gradient.(Ref(itp), idx) ==
        collect(Interpolations.gradient.(Ref(jlitp), idx)) ==
        collect(Interpolations.gradient.(Ref(jlitp), jlidx))

    sitp = scale(itp, A_x)
    jlsitp = jl(sitp)
    idx = 1.0:0.4:39.0
    jlidx = jl(collect(idx))
    @test sitp.(idx) == collect(jlsitp.(idx)) == collect(jlsitp.(jlidx))
    @test Interpolations.gradient.(Ref(sitp), idx) ==
        collect(Interpolations.gradient.(Ref(jlsitp), idx)) ==
        collect(Interpolations.gradient.(Ref(jlsitp), jlidx))


    esitp = extrapolate(sitp, Flat())
    jlesitp = jl(esitp)
    idx = -1.0:0.84:41.0
    jlidx = jl(collect(idx))
    @test esitp.(idx) == collect(jlesitp.(idx)) == collect(jlesitp.(jlidx))
    @test Interpolations.gradient.(Ref(esitp), idx) ==
        collect(Interpolations.gradient.(Ref(jlesitp), idx)) ==
        collect(Interpolations.gradient.(Ref(jlesitp), jlidx))

    esitp = extrapolate(sitp, 0.0)
    jlesitp = jl(esitp)
    idx = -1.0:0.84:41.0
    jlidx = jl(collect(idx))
    @test esitp.(idx) == collect(jlesitp.(idx)) == collect(jlesitp.(jlidx))
    @test Interpolations.gradient.(Ref(esitp), idx) ==
        collect(Interpolations.gradient.(Ref(jlesitp), idx)) ==
        collect(Interpolations.gradient.(Ref(jlesitp), jlidx))
end

@testset "2d GPU Interpolation" begin
    A_x = 1.0:2.0:40.0
    A = [log(x + y) for x in A_x, y in 1.0:2.0:40.0]
    itp = interpolate(A, (BSpline(Cubic(Line(OnGrid()))), BSpline(Linear())))
    jlitp = jl(itp)
    idx = 2.0:0.17:19.0
    jlidx = jl(collect(idx))
    @test itp.(idx, idx') == collect(jlitp.(idx, idx')) == collect(jlitp.(jlidx, jlidx'))
    @test Interpolations.gradient.(Ref(itp), idx, idx') ==
        collect(Interpolations.gradient.(Ref(jlitp), idx, idx')) ==
        collect(Interpolations.gradient.(Ref(jlitp), jlidx, jlidx'))
    @test Interpolations.hessian.(Ref(itp), idx, idx') ==
        collect(Interpolations.hessian.(Ref(jlitp), idx, idx')) ==
        collect(Interpolations.hessian.(Ref(jlitp), jlidx, jlidx'))

    sitp = scale(itp, A_x, A_x)
    jlsitp = jl(sitp)
    idx = 1.0:0.4:39.0
    jlidx = jl(collect(idx))
    @test sitp.(idx, idx') == collect(jlsitp.(idx, idx')) == collect(jlsitp.(jlidx, jlidx'))
    @test Interpolations.gradient.(Ref(sitp), idx, idx') ==
        collect(Interpolations.gradient.(Ref(jlsitp), idx, idx')) ==
        collect(Interpolations.gradient.(Ref(jlsitp), jlidx, jlidx'))
    @test Interpolations.hessian.(Ref(sitp), idx, idx') ==
        collect(Interpolations.hessian.(Ref(jlsitp), idx, idx')) ==
        collect(Interpolations.hessian.(Ref(jlsitp), jlidx, jlidx'))

    esitp = extrapolate(sitp, Flat())
    jlesitp = jl(esitp)
    idx = -1.0:0.84:41.0
    jlidx = jl(collect(idx))
    @test esitp.(idx, idx') == collect(jlesitp.(idx, idx')) == collect(jlesitp.(jlidx, jlidx'))
    # Interpolations.gradient for `extrapolation` is currently broken under CUDA
    @test Interpolations.gradient.(Ref(esitp), idx, idx') ==
        collect(Interpolations.gradient.(Ref(jlesitp), idx, idx')) ==
        collect(Interpolations.gradient.(Ref(jlesitp), jlidx, jlidx'))

    esitp = extrapolate(sitp, 0.0)
    jlesitp = jl(esitp)
    idx = -1.0:0.84:41.0
    jlidx = jl(collect(idx))
    @test esitp.(idx, idx') == collect(jlesitp.(idx, idx')) == collect(jlesitp.(jlidx, jlidx'))
    # Interpolations.gradient for `extrapolation` is currently broken under CUDA
    @test Interpolations.gradient.(Ref(esitp), idx, idx') ==
        collect(Interpolations.gradient.(Ref(jlesitp), idx, idx')) ==
        collect(Interpolations.gradient.(Ref(jlesitp), jlidx, jlidx'))
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
    itp1 = interpolate(Vector.((-1.0:0.02:1.0, -1.0:0.02:1.0)), randn(101, 101), Gridded(Linear()))
    itp2 = interpolate((-1.0:0.02:1.0, -1.0:0.02:1.0), randn(101, 101), Gridded(Linear()))
    idx = -1.0:0.01:1.0
    jlidx = jl(collect(idx))
    @test itp1.(idx, idx') == collect(jl(itp1).(idx, idx')) == collect(jl(itp1).(jlidx, jlidx'))
    @test itp2.(idx, idx') == collect(jl(itp2).(idx, idx')) == collect(jl(itp2).(jlidx, jlidx'))
end

@testset "eltype after adaption" begin
    A_x = 1.0:2.0:40.0
    A = [log(x) for x in A_x]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    @test eltype(adapt(Array{Float32}, itp)) === Float32
    @test eltype(adapt(Array{Real}, itp)) === Float64
    @test eltype(adapt(Array{Float32}, scale(itp, A_x))) === Float32
    @test eltype(adapt(Array{Float32}, extrapolate(scale(itp, A_x), Flat()))) === Float32
    @test eltype(adapt(Array{Float32}, extrapolate(scale(itp, A_x), 0.0))) === Float32
    itp = interpolate((-1:0.2:1, -1:0.2:1), randn(11, 11), Gridded(Linear()))
    @test eltype(adapt(Array{Float32}, itp)) === Float32
    itp = interpolate((1.0:0.0, 1.:0.), randn(0, 0), Gridded(Linear()))
    @test eltype(adapt(Array{Real}, itp)) isa DataType # we don't care the result
end
