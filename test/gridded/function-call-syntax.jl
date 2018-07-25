module GriddedFunctionCallSyntax

using Interpolations, Compat.Test
using Compat: range

for D in (Constant, Linear), G in (OnCell, OnGrid)
    ## 1D
    a = rand(5)
    knots = (collect(range(1, stop=length(a), length=length(a))),)
    itp = @inferred(interpolate(knots, a, Gridded(D())))
    @inferred(getindex(itp, 2))
    @inferred(getindex(itp, CartesianIndex((2,))))
    for i = 2:length(a)-1
        @test itp(i) ≈ a[i]
        @test itp(CartesianIndex((i,))) ≈ a[i]
    end
    @inferred(getindex(itp, knots...))
    @test itp[knots...] ≈ a
    # compare scalar indexing and vector indexing
    x = knots[1] .+ 0.1
    v = itp(x)
    for i = 1:length(x)
        @test v[i] ≈ itp(x[i])
    end
    # check the fallback vector indexing
    x = [2.3,2.2]   # non-increasing order
    v = itp[x]
    for i = 1:length(x)
        @test v[i] ≈ itp(x[i])
    end
    # compare against BSpline
    itpb = @inferred(interpolate(a, BSpline(D()), G()))
    for x in range(1.1, stop=4.9, length=101)
        @test itp(x) ≈ itpb(x)
    end

    ## 2D
    A = rand(6,5)
    knots = (collect(range(1, stop=size(A,1), length=size(A,1))),collect(range(1, stop=size(A,2), length=size(A,2))))
    itp = @inferred(interpolate(knots, A, Gridded(D())))
    @test parent(itp) === A
    @inferred(getindex(itp, 2, 2))
    @inferred(getindex(itp, CartesianIndex((2,2))))
    for j = 2:size(A,2)-1, i = 2:size(A,1)-1
        @test itp(i,j) ≈ A[i,j]
        @test itp(CartesianIndex((i,j))) ≈ A[i,j]
    end
    @test itp[knots...] ≈ A
    @inferred(getindex(itp, knots...))
    # compare scalar indexing and vector indexing
    x, y = knots[1] .+ 0.1, knots[2] .+ 0.6
    v = itp(x,y)
    for j = 1:length(y), i = 1:length(x)
        @test v[i,j] ≈ itp(x[i],y[j])
    end
    # check the fallback vector indexing
    x = [2.3,2.2]   # non-increasing order
    y = [3.5,2.8]
    v = itp[x,y]
    for j = 1:length(y), i = 1:length(x)
        @test v[i,j] ≈ itp(x[i],y[j])
    end
    # compare against BSpline
    itpb = @inferred(interpolate(A, BSpline(D()), G()))
    for y in range(1.1, stop=5.9, length=101), x in range(1.1, stop=4.9, length=101)
        @test itp(x,y) ≈ itpb(x,y)
    end

    A = rand(8,20)
    knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
    itp = interpolate(knots, A, Gridded(D()))
    @test itp(4,1.2) ≈ A[2,6]
end

end
