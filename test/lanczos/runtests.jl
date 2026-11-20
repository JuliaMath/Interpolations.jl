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

@testset "Lanczos4OpenCVFaithful" begin
    # Strict input element-type enforcement (only UInt8 and Float32) & type stability
    @test (@inferred Lanczos4OpenCVFaithful()) === Lanczos4OpenCVFaithful{false}()
    @test Lanczos4OpenCVFaithful(simd=true) === Lanczos4OpenCVFaithful{true}()
    @test Lanczos4OpenCVFaithful(simd=false) === Lanczos4OpenCVFaithful{false}()
    @test Interpolations.value_weights(Lanczos4OpenCVFaithful(), 0.0f0)[4] === 1.0f0
    @test Interpolations.value_weights(Lanczos4OpenCV(), 0.0)[4] == 1.0
    @test_throws ArgumentError Lanczos4OpenCVFaithful{Int}()
    @test_throws ArgumentError interpolate(zeros(Float64, 8, 8), Lanczos4OpenCVFaithful())
    @test_throws ArgumentError interpolate(zeros(Float16, 8, 8), Lanczos4OpenCVFaithful())
    @test_throws ArgumentError interpolate(zeros(Int, 8, 8), Lanczos4OpenCVFaithful())
    @test_throws ArgumentError interpolate(zeros(UInt16, 8, 8), Lanczos4OpenCVFaithful())

    # Verify LUT builder matches precompiled const tables
    t1d, t1d_i, t2d_f, t2d_i = Interpolations._build_lanczos4_ocv_tables()
    @test t1d == Interpolations._LANCZOS4_OCV_TAB1D
    @test t1d_i == Interpolations._LANCZOS4_OCV_TAB1D_I
    @test t2d_f == Interpolations._LANCZOS4_OCV_TAB_F
    @test t2d_i == Interpolations._LANCZOS4_OCV_TAB_I

    for T in (UInt8, Float32)
        v = T.(collect(1:100))
        X1 = [v; reverse(v)[2:end]]
        itp1 = interpolate(X1, Lanczos4OpenCVFaithful())
        itp1_simd = interpolate(X1, Lanczos4OpenCVFaithful(simd=true))
        @test eltype(itp1) === T
        @test eltype(coefficients(itp1)) === T
        @test getknots(itp1) == axes(itp1) == axes(X1)
        @test coefficients(itp1)[axes(itp1)...] == X1
        @test itpflag(itp1) isa Lanczos4OpenCVFaithful
        @test degree(itpflag(itp1)) == 4
        @test size(itp1) == size(X1)
        @test lbounds(itp1) == (firstindex(X1),)
        @test ubounds(itp1) == (lastindex(X1),)
        @test [itp1(i) for i in eachindex(X1)] == X1
        @test [itp1_simd(i) for i in eachindex(X1)] == X1
        @test isapprox(itp1_simd(10.25f0), itp1(10.25f0); atol=1e-3)
        @test itp1(10.25f0) === itp1(10.25)
    end

    # 3D fallback path (positions / value_weights / N >= 3 UInt8 rounding)
    a3_u8 = reshape(UInt8[(i * 17) % 256 for i in 1:512], 8, 8, 8)
    itp3_u8 = interpolate(a3_u8, Lanczos4OpenCVFaithful())
    itp3_f32 = interpolate(Float32.(a3_u8), Lanczos4OpenCVFaithful())
    @test itp3_u8(4, 4, 4) === a3_u8[4, 4, 4]
    @test itp3_u8(4.25f0, 4.5f0, 4.125) isa UInt8
    @test isapprox(Float32(itp3_u8(4.25f0, 4.5f0, 4.125)), itp3_f32(4.25f0, 4.5f0, 4.125); atol=1.0f0)

    # 2D exact knot recovery & symmetry for both UInt8 and Float32
    p_u8 = UInt8[(c * 73 + r * 37 + r * c * 11) % 256 for c in 0:11, r in 0:11]
    p_f32 = Float32[sin(0.7 * r + 0.3 * c) * 100.0 for c in 0:11, r in 0:11]
    itp_u8 = interpolate(p_u8, Lanczos4OpenCVFaithful())
    itp_f32 = interpolate(p_f32, Lanczos4OpenCVFaithful())
    @test [itp_u8(i, j) for i in axes(p_u8, 1), j in axes(p_u8, 2)] == p_u8
    @test [itp_f32(i, j) for i in axes(p_f32, 1), j in axes(p_f32, 2)] == p_f32

    # 1/32-pixel quantization check: offsets within 1/64 of each other in the same bin agree identically
    @test itp_f32(5.25f0, 6.5f0) === itp_f32(5.25f0 + 0.001f0, 6.5f0 - 0.001f0)
    @test itp_u8(5.25f0, 6.5f0) === itp_u8(5.25f0 + 0.001f0, 6.5f0 - 0.001f0)

    # Stencil shift at >= 31.5/32 rounds up to next integer knot
    @test itp_f32(4.0 + 31.9 / 32.0, 5.0 + 31.9 / 32.0) === p_f32[5, 6]
    @test itp_u8(4.0 + 31.9 / 32.0, 5.0 + 31.9 / 32.0) === p_u8[5, 6]

    # Exact bit-for-bit parity against OpenCV 4.11 cv::remap(..., INTER_LANCZOS4, BORDER_CONSTANT, 0)
    qx = (0.0f0, 3.25f0, 5.5f0 + 0.5f0 / 32.0f0, 7.0f0 + 31.9f0 / 32.0f0, 11.0f0)
    qy = (0.0f0, 4.125f0, 6.75f0, 2.0f0 + 15.5f0 / 32.0f0, 10.5f0)
    @test [itp_u8(Float64(x) + 1.0, Float64(y) + 1.0) for (x, y) in zip(qx, qy)] ==
          UInt8[0, 201, 210, 131, 166]
    @test [itp_f32(Float64(x) + 1.0, Float64(y) + 1.0) for (x, y) in zip(qx, qy)] ==
          Float32[0.0, -0x1.08eaacp6, 0x1.237146p3, -0x1.4bea2p6, -0x1.a2ca0ap6]

    # Arbitrarily large (> 32k SHRT_MAX and > 524k Float32 mantissa) coordinate / axis support
    big_off = 2_000_000
    oa_u8 = OffsetArrays.OffsetArray(p_u8, (1 + big_off):(12 + big_off), (1 + big_off):(12 + big_off))
    oa_f32 = OffsetArrays.OffsetArray(p_f32, (1 + big_off):(12 + big_off), (1 + big_off):(12 + big_off))
    itp_big_u8 = interpolate(oa_u8, Lanczos4OpenCVFaithful())
    itp_big_f32 = interpolate(oa_f32, Lanczos4OpenCVFaithful())
    @test [itp_big_u8(Float64(x) + 1.0 + big_off, Float64(y) + 1.0 + big_off) for (x, y) in zip(qx, qy)] ==
          UInt8[0, 201, 210, 131, 166]
    @test [itp_big_f32(Float64(x) + 1.0 + big_off, Float64(y) + 1.0 + big_off) for (x, y) in zip(qx, qy)] ==
          Float32[0.0, -0x1.08eaacp6, 0x1.237146p3, -0x1.4bea2p6, -0x1.a2ca0ap6]
    # Also test a > 32,768-length 1D array directly (where OpenCV cv::remap fails SHRT_MAX)
    big_vec = UInt8[(i * 13) % 256 for i in 1:40_000]
    itp_40k = interpolate(big_vec, Lanczos4OpenCVFaithful())
    @test itp_40k(39_999) === big_vec[39_999]
    @test itp_40k(35_000.25) === interpolate(big_vec[34_995:35_006], Lanczos4OpenCVFaithful())(6.25)

    # Relaxed SIMD FMA mode (simd=true): exact at knots, within a couple ULPs of simd=false on Float32
    itp_f32_simd = @inferred interpolate(p_f32, Lanczos4OpenCVFaithful(simd=true))
    @test [itp_f32_simd(i, j) for i in axes(p_f32, 1), j in axes(p_f32, 2)] == p_f32
    for (x, y) in zip(qx, qy)
        v_exact = @inferred itp_f32(Float64(x) + 1.0, Float64(y) + 1.0)
        v_simd = @inferred itp_f32_simd(Float64(x) + 1.0, Float64(y) + 1.0)
        @test isapprox(v_simd, v_exact; atol=4 * eps(Float32(max(1.0f0, abs(v_exact)))))
    end

    # Type stability & zero-allocation checks across direct, scaled, and filled-extrapolation calls
    _alloc2(f::F, x, y) where {F} = (f(x, y); @allocated f(x, y))
    for (itp2, z) in ((itp_u8, 0x00), (itp_f32, 0.0f0), (itp_f32_simd, 0.0f0))
        T = eltype(itp2)
        etp2 = @inferred extrapolate(itp2, z)
        sitp2 = @inferred scale(itp2, 1:12, 1:12)
        estp2 = @inferred extrapolate(sitp2, z)
        for (x, y) in ((5, 6), (5.25f0, 6.75f0), (5.25, 6.75))
            @test (@inferred itp2(x, y)) isa T
            @test (@inferred etp2(x, y)) isa T
            @test (@inferred etp2(-3.0, 99.0)) === z
            @test (@inferred sitp2(x, y)) isa T
            @test (@inferred estp2(x, y)) isa T
            @test _alloc2(itp2, x, y) == 0
            @test _alloc2(etp2, x, y) == 0
            @test _alloc2(sitp2, x, y) == 0
            @test _alloc2(estp2, x, y) == 0
        end
    end
end

end
