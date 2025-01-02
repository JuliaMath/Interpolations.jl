using Test

function convolution_test()
    N = 3

    ### 1D
    range_1d = range(0.0, stop=1.0, length=N)

    # random uniformly spaced 1D data
    vals_1d_rand = rand(N)
    itp_1d_rand = convolution_interpolation(range_1d, vals_1d_rand, extrapolation_bc=Interpolations.Line())
    @test itp_1d_rand.(range_1d) ≈ vals_1d_rand
    
    # match mean of linear uniformly spaced 1D data
    vals_1d_linear = range(0.0, stop=1.0, length=N)
    vals_1d_linear_mean =[(vals_1d_linear[i]+vals_1d_linear[i+1])/2 for i in 1:N-1]
    itp_1d_linear = convolution_interpolation(range_1d, vals_1d_linear, extrapolation_bc=Interpolations.Line())
    @test itp_1d_linear.(range(vals_1d_linear_mean[1], stop=vals_1d_linear_mean[end], length=N-1)) ≈ vals_1d_linear_mean

    # extrapolation
    etp_1d_linear = convolution_interpolation(range_1d, vals_1d_linear, extrapolation_bc=Interpolations.Line())
    @test etp_1d_linear.([-0.2, -0.1, 1.1, 1.2]) ≈ [-0.2, -0.1, 1.1, 1.2]
    etp_1d_flat = convolution_interpolation(range_1d, vals_1d_linear, extrapolation_bc=Interpolations.Flat())
    @test etp_1d_flat.([-0.2, -0.1, 1.1, 1.2]) ≈ [0.0, 0.0, 1.0, 1.0]

    ### 2D
    range_2d = range(0.0, stop=1.0, length=N)
    
    # random uniformly spaced 2D data
    vals_2d_rand = rand(N,N)
    itp_2d_rand = convolution_interpolation((range_1d, range_2d), vals_2d_rand, extrapolation_bc=Interpolations.Line())
    @test [itp_2d_rand(range_1d[i], range_2d[j]) for i in 1:N for j in 1:N] ≈ [vals_2d_rand[i,j] for i in 1:N for j in 1:N]

    # mean of linear uniformly spaced 2D data
    f_2d_linear(i,j) = (i-1)/(N-1)+(j-1)/(N-1)
    vals_2d_linear = Matrix{Float64}(undef, N, N)
    [vals_2d_linear[i,j] = f_2d_linear(i,j) for i in 1:N for j in 1:N]
    vals_2d_linear_mean =[(f_2d_linear(i,j)+f_2d_linear(i+1,j)+f_2d_linear(i,j+1)+f_2d_linear(i+1,j+1))/4 for i in 1:N-1 for j in 1:N-1]
    itp_2d_linear = convolution_interpolation((range_1d, range_2d,), vals_2d_linear, extrapolation_bc=Interpolations.Line())
    @test [itp_2d_linear((i-1/2)/(N-1),(j-1/2)/(N-1)) for i in 1:N-1 for j in 1:N-1][:] ≈ vals_2d_linear_mean

    # extrapolation
    etp_2d_linear = convolution_interpolation((range_1d, range_2d,), vals_2d_linear, extrapolation_bc=Interpolations.Line())
    @test [etp_2d_linear(i, i) for i in [-0.2, -0.1, 1.1, 1.2]] ≈ [-0.4, -0.2, 2.2, 2.4]
    etp_2d_flat = convolution_interpolation((range_1d, range_2d,), vals_2d_linear, extrapolation_bc=Interpolations.Flat())
    @test [etp_2d_flat(i, i) for i in [-0.2, -0.1, 1.1, 1.2]] ≈ [0.0, 0.0, 2.0, 2.0]
    
    ### 3D
    range_3d = range(0.0, stop=1.0, length=N)

    # random uniformly spaced 3D data
    vals_3d_rand = rand(N,N,N)
    itp_3d_rand = convolution_interpolation((range_1d, range_2d, range_3d), vals_3d_rand, extrapolation_bc=Interpolations.Line())
    @test [itp_3d_rand(range_1d[i], range_2d[j], range_3d[k]) for i in 1:N for j in 1:N for k in 1:N] ≈ [vals_3d_rand[i,j,k] for i in 1:N for j in 1:N for k in 1:N]

    # mean of linear uniformly spaced 3D data
    f_3d_linear(i,j,k) = (i-1)/(N-1)+(j-1)/(N-1)+(k-1)/(N-1)
    vals_3d_linear = Array{Float64}(undef, N, N, N)
    [vals_3d_linear[i,j,k] = f_3d_linear(i,j,k) for i in 1:N for j in 1:N for k in 1:N]
    vals_3d_linear_mean =[(f_3d_linear(i,j,k)+f_3d_linear(i+1,j,k)+f_3d_linear(i,j+1,k)+f_3d_linear(i,j,k+1)+
                            f_3d_linear(i+1,j+1,k+1)+f_3d_linear(i,j+1,k+1)+f_3d_linear(i+1,j,k+1)+f_3d_linear(i+1,j+1,k))/8 for i in 1:N-1 for j in 1:N-1 for k in 1:N-1]
    itp_3d_linear = convolution_interpolation((range_1d, range_2d, range_3d), vals_3d_linear, extrapolation_bc=Interpolations.Line())
    @test [itp_3d_linear((i-1/2)/(N-1),(j-1/2)/(N-1),(k-1/2)/(N-1)) for i in 1:N-1 for j in 1:N-1 for k in 1:N-1][:] ≈ vals_3d_linear_mean

    # extrapolation
    etp_3d_linear = convolution_interpolation((range_1d, range_2d, range_3d), vals_3d_linear, extrapolation_bc=Interpolations.Line())
    @test isapprox([etp_3d_linear(i, i, i) for i in [-0.2, -0.1, 1.1, 1.2]], [-0.6, -0.3, 3.3, 3.6], atol=1e-3)
    etp_3d_flat = convolution_interpolation((range_1d, range_2d, range_3d), vals_3d_linear, extrapolation_bc=Interpolations.Flat())
    @test [etp_3d_flat(i, i, i) for i in [-0.2, -0.1, 1.1, 1.2]] ≈ [0.0, 0.0, 3.0, 3.0]

    ### 4D
    range_4d = range(0.0, stop=1.0, length=N)

    # random uniformly spaced 4D data
    vals_4d_rand = rand(N,N,N,N)
    itp_4d_rand = convolution_interpolation((range_1d, range_2d, range_3d, range_4d), vals_4d_rand, extrapolation_bc=Interpolations.Line())
    @test [itp_4d_rand(range_1d[i], range_2d[j], range_3d[k], range_4d[l]) for i in 1:N for j in 1:N for k in 1:N for l in 1:N] ≈ [vals_4d_rand[i,j,k,l] for i in 1:N for j in 1:N for k in 1:N for l in 1:N]

    # mean of linear uniformly spaced 4D data
    f_4d_linear(i,j,k,l) = (i-1)/(N-1)+(j-1)/(N-1)+(k-1)/(N-1)+(l-1)/(N-1)
    vals_4d_linear = Array{Float64}(undef, N, N, N, N)
    [vals_4d_linear[i,j,k,l] = f_4d_linear(i,j,k,l) for i in 1:N for j in 1:N for k in 1:N for l in 1:N]
    function hypercube_mean(f_4d_linear, N)
        [(
            f_4d_linear(i, j, k, l) +
            f_4d_linear(i+1, j, k, l) +
            f_4d_linear(i, j+1, k, l) +
            f_4d_linear(i, j, k+1, l) +
            f_4d_linear(i, j, k, l+1) +
            f_4d_linear(i+1, j+1, k, l) +
            f_4d_linear(i+1, j, k+1, l) +
            f_4d_linear(i+1, j, k, l+1) +
            f_4d_linear(i, j+1, k+1, l) +
            f_4d_linear(i, j+1, k, l+1) +
            f_4d_linear(i, j, k+1, l+1) +
            f_4d_linear(i+1, j+1, k+1, l) +
            f_4d_linear(i+1, j+1, k, l+1) +
            f_4d_linear(i+1, j, k+1, l+1) +
            f_4d_linear(i, j+1, k+1, l+1) +
            f_4d_linear(i+1, j+1, k+1, l+1)
        ) / 16 for i in 1:N-1, j in 1:N-1, k in 1:N-1, l in 1:N-1]
    end
    vals_4d_linear_mean = hypercube_mean(f_4d_linear, N)
    itp_4d_linear = convolution_interpolation((range_1d, range_2d, range_3d, range_4d), vals_4d_linear, extrapolation_bc=Interpolations.Line())
    @test [itp_4d_linear((i-1/2)/(N-1),(j-1/2)/(N-1),(k-1/2)/(N-1),(l-1/2)/(N-1)) for i in 1:N-1 for j in 1:N-1 for k in 1:N-1 for l in 1:N-1] ≈ vals_4d_linear_mean[:]

    # test extrapolation bc
    etp_4d_linear = convolution_interpolation((range_1d, range_2d, range_3d, range_4d), vals_4d_linear, extrapolation_bc=Interpolations.Line())
    @test [etp_4d_linear( i, i, i, i) for i in [-0.2, -0.1, 1.1, 1.2]] ≈ [-0.8, -0.4, 4.4, 4.8]
    etp_4d_flat = convolution_interpolation((range_1d, range_2d, range_3d, range_4d), vals_4d_linear, extrapolation_bc=Interpolations.Flat())
    @test [etp_4d_flat( i, i, i, i) for i in [-0.2, -0.1, 1.1, 1.2]] ≈ [0.0, 0.0, 4.0, 4.0]

end