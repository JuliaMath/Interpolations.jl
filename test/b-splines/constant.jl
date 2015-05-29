module ConstantTests

using Interpolations
using Base.Test

# Instantiation
N1 = 10
A1 = rand(Float64, N1) * 100
A2 = rand(Float64, N1, N1) * 100
A3 = rand(Float64, N1, N1, N1) * 100

for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
    itp1c = @inferred(constructor(copier(A1), BSpline(Constant), OnCell))
    itp1g = @inferred(constructor(copier(A1), BSpline(Constant), OnGrid))
    itp2c = @inferred(constructor(copier(A2), BSpline(Constant), OnCell))
    itp2g = @inferred(constructor(copier(A2), BSpline(Constant), OnGrid))
    itp3c = @inferred(constructor(copier(A3), BSpline(Constant), OnCell))
    itp3g = @inferred(constructor(copier(A3), BSpline(Constant), OnGrid))

    # Evaluation on provided data points
    # 1D
    for i in 1:length(A1)
        @test A1[i] == itp1c[i] == itp1g[i]
        @test A1[i] == itp1c[convert(Float64,i)] == itp1g[convert(Float64,i)]
    end
    @test @inferred(size(itp1c)) == size(A1)
    @test @inferred(size(itp1g)) == size(A1)
    # 2D
    for i in 1:N1, j in 1:N1
        @test A2[i,j] == itp2c[i,j] == itp2g[i,j]
        @test A2[i,j] == itp2c[convert(Float64,i),convert(Float64,j)] == itp2g[convert(Float64,i),convert(Float64,j)]
    end
    @test @inferred(size(itp2c)) == size(A2)
    @test @inferred(size(itp2g)) == size(A2)
    # 3D
    for i in 1:N1, j in 1:N1, k in 1:N1
        @test A3[i,j,k] == itp3c[i,j,k] == itp3g[i,j,k]
        @test A3[i,j,k] == itp3c[convert(Float64,i),convert(Float64,j),convert(Float64,k)] == itp3g[convert(Float64,i),convert(Float64,j),convert(Float64,k)]
    end
    @test @inferred(size(itp3c)) == size(A3)
    @test @inferred(size(itp3g)) == size(A3)

    # Evaluation between data points
    for i in 2:N1-1
        @test A1[i] == itp1c[i+.3] == itp1g[i+.3] == itp1c[i-.3] == itp1g[i-.3]
    end
    # 2D
    for i in 2:N1-1, j in 2:N1-1
        @test A2[i,j] == itp2c[i+.4,j-.3] == itp2g[i+.4,j-.3]
    end
    # 3D
    for i in 2:N1-1, j in 2:N1-1, k in 2:N1-1
        @test A3[i,j,k] == itp3c[i+.4,j-.3,k+.1] == itp3g[i+.4,j-.3,k+.2]
    end

    # Edge behavior
    @test A1[1] == itp1c[.7]
    @test A1[N1] == itp1c[N1+.3]
end

end