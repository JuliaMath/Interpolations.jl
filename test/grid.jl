module GridTests

using Interpolations, Base.Test

# On-grid values
A = randn(4,10)
const EPS = sqrt(eps())
for it in (Constant(OnCell()), Linear(OnGrid()), Quadratic(Free(),OnCell()))
    for eb in (ExtrapNaN(), ExtrapError())
        itp = Interpolation(A, it, eb)
        for i = 1:size(A,1)
            for j = 1:size(A,2)
                @test_approx_eq_eps itp[i,j] A[i,j] EPS
            end
        end
        # v = itp[1:size(A,1), 1:size(A,2)]
        # @assert all(abs(v - A) .< EPS)
    end
end

A = randn(4,5,4)
for it in (Constant(OnCell()), Linear(OnGrid()), Quadratic(Free(),OnCell()))
    for eb in (ExtrapNaN(), ExtrapError())
        itp = Interpolation(A, it, eb)
        for k = 1:size(A,3), j = 1:size(A,2), i = 1:size(A,1)
            @test_approx_eq_eps itp[i,j,k] A[i,j,k] EPS
        end
        # v = itp[1:size(A,1), 1:size(A,2), 1:size(A,3)]
        # @assert all(abs(v - A) .< EPS)
    end
end
A = randn(4,5,4,3)
for it in (Constant(OnCell()), Linear(OnGrid()), Quadratic(Free(),OnCell()))
    for eb in (ExtrapNaN(), ExtrapError())
        itp = Interpolation(A, it, eb)
        for i = 1:size(A,1), j = 1:size(A,2), k = 1:size(A,3), l = 1:size(A,4)
            @test_approx_eq_eps itp[i,j,k,l] A[i,j,k,l] EPS
        end
        # v = itp[1:size(A,1), 1:size(A,2), 1:size(A,3), 1:size(A,4)]
        # @assert all(abs(v - A) .< EPS)
    end
end

A = float([1:4])
for it in (Constant(OnCell()), Linear(OnGrid()), Quadratic(Free(),OnCell()))
    itp = Interpolation(A, it, ExtrapError())
    @test_throws BoundsError itp[-0.8]
end
for it in (Constant(OnCell()), Linear(OnGrid()), Quadratic(Free(),OnCell()))
    itp = Interpolation(A, it, ExtrapNaN())
    @test isnan(itp[-0.8])
end

end
