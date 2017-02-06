module MixedTests

using Interpolations, Base.Test

N = 10

for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
    A2 = rand(Float64, N, N) * 100
    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        itp_a = @inferred(constructor(copier(A2), (BSpline(Linear()), BSpline(Quadratic(BC()))), GT()))
        itp_b = @inferred(constructor(copier(A2), (BSpline(Quadratic(BC())), BSpline(Linear())), GT()))
        @test @inferred(size(itp_a)) == size(A2)
        @test @inferred(size(itp_b)) == size(A2)

        for j = 2:N-1, i = 2:N-1
            @test_approx_eq_eps itp_a[i,j] A2[i,j] sqrt(eps(A2[i,j]))
            @test_approx_eq_eps itp_b[i,j] A2[i,j] sqrt(eps(A2[i,j]))
        end

        for i = 1:10
            dx, dy = rand(), rand()
            @test_approx_eq itp_a[2+dx,2] (1-dx)*A2[2,2]+dx*A2[3,2]
            @test_approx_eq itp_b[2,2+dy] (1-dy)*A2[2,2]+dy*A2[2,3]
        end
    end
end

# AbstractArrays
function copyshared(A)
    B = SharedArray{eltype(A)}( size(A))
    copy!(B, A)
end

for (constructor, copier) in ((interpolate, x->x), (interpolate!, copyshared))
    A2 = SharedArray{Float64}( (N,N), init=A->rand!(A))
    for i = 1:length(A2)
        A2[i] *= 100
    end
    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        itp_a = @inferred(constructor(copier(A2), (BSpline(Linear()), BSpline(Quadratic(BC()))), GT()))
        itp_b = @inferred(constructor(copier(A2), (BSpline(Quadratic(BC())), BSpline(Linear())), GT()))
        if constructor == interpolate!
            @test isa(itp_a.coefs, SharedArray)
            @test isa(itp_b.coefs, SharedArray)
        end
        @test @inferred(size(itp_a)) == size(A2)
        @test @inferred(size(itp_b)) == size(A2)

        for j = 2:N-1, i = 2:N-1
            @test_approx_eq_eps itp_a[i,j] A2[i,j] sqrt(eps(A2[i,j]))
            @test_approx_eq_eps itp_b[i,j] A2[i,j] sqrt(eps(A2[i,j]))
        end

        for i = 1:10
            dx, dy = rand(), rand()
            @test_approx_eq itp_a[2+dx,2] (1-dx)*A2[2,2]+dx*A2[3,2]
            @test_approx_eq itp_b[2,2+dy] (1-dy)*A2[2,2]+dy*A2[2,3]
        end
    end
end

end