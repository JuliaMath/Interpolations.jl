module ConstantTests

using Interpolations
using Test

# Instantiation
N1 = 10
A1 = rand(Float64, N1) * 100
A2 = rand(Float64, N1, N1) * 100
A3 = rand(Float64, N1, N1, N1) * 100

getindexib(itp, i...) = @inbounds itp[i...]
callib(itp, i...) = @inbounds itp(i...)

for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
    itp1c = @inferred(constructor(copier(A1), BSpline(Constant()), OnCell()))
    itp1g = @inferred(constructor(copier(A1), BSpline(Constant()), OnGrid()))
    itp2c = @inferred(constructor(copier(A2), BSpline(Constant()), OnCell()))
    itp2g = @inferred(constructor(copier(A2), BSpline(Constant()), OnGrid()))
    itp3c = @inferred(constructor(copier(A3), BSpline(Constant()), OnCell()))
    itp3g = @inferred(constructor(copier(A3), BSpline(Constant()), OnGrid()))

    @test parent(itp1c) === itp1c.coefs

    # Evaluation on provided data points
    for (itp, A) in ((itp1c, A1), (itp1g, A1), (itp2c, A2), (itp2g, A2), (itp3c, A3), (itp3g, A3))
        for i in eachindex(itp)
            @test A[i] == itp[i] == itp[Tuple(i)...] == itp(i) == itp(float.(Tuple(i))...)
        end
        @test @inferred(axes(itp)) == axes(A)
        @test @inferred(size(itp)) == size(A)
        # @test @inferred(axes(itp)) == (contructor == interpolate ? (1:N1) : (2:N1-1))
        # @test @inferred(size(itp)) == (contructor == interpolate ? N1 : N1-2)
    end
    for itp in (itp1c, itp1g)
        for i = 1:N1
            @test itp[i,1] == A1[i]   # used in the AbstractArray display infrastructure
            @test_throws BoundsError itp[i,2]
            @test_broken getindexib(itp, i, 2) == A1[i]
            @test_broken callib(itp, i, 2) == A1[i]
        end
    end

    # Evaluation between data points (tests constancy)
    for i in 2:N1-1
        @test A1[i] == itp1c(i+.3) == itp1g(i+.3) == itp1c(i-.3) == itp1g(i-.3)
    end
    # 2D
    for i in 2:N1-1, j in 2:N1-1
        @test A2[i,j] == itp2c(i+.4,j-.3) == itp2g(i+.4,j-.3)
    end
    # 3D
    for i in 2:N1-1, j in 2:N1-1, k in 2:N1-1
        @test A3[i,j,k] == itp3c(i+.4,j-.3,k+.1) == itp3g(i+.4,j-.3,k+.2)
    end

    # Edge behavior
    @test A1[1] == itp1c(.7)
    @test A1[N1] == itp1c(N1+.3)
end

end
