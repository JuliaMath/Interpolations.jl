@testset "Mixed" begin
    N = 10

    for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
        A2 = rand(Float64, N, N) * 100
        for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
            itp_a = @inferred(constructor(copier(A2), (BSpline(Linear()), BSpline(Quadratic(BC(GT()))))))
            itp_b = @inferred(constructor(copier(A2), (BSpline(Quadratic(BC(GT()))), BSpline(Linear()))))
            isfullsize = constructor == interpolate || BC==Periodic
            if isfullsize
                @test @inferred(size(itp_a)) == size(A2)
                @test @inferred(size(itp_b)) == size(A2)
                @test @inferred(axes(itp_a)) == axes(A2)
                @test @inferred(axes(itp_b)) == axes(A2)
            else
                @test @inferred(size(itp_a)) == (N, N-2)
                @test @inferred(size(itp_b)) == (N-2, N)
                @test @inferred(axes(itp_a)) == (1:N, 2:N-1)
                @test @inferred(axes(itp_b)) == (2:N-1, 1:N)
            end
            @test_throws ArgumentError parent(itp_a)
            @test_throws ArgumentError parent(itp_b)

            for i in eachindex(itp_a)
                @test itp_a[i] ≈ A2[i] atol=sqrt(eps(A2[i]))
            end
            for i in eachindex(itp_b)
                @test itp_b[i] ≈ A2[i] atol=sqrt(eps(A2[i]))
            end

            for i = 1:10
                dx, dy = rand(), rand()
                @test itp_a(2 + dx,2) ≈ (1 - dx) * A2[2,2] + dx * A2[3,2]
                @test itp_b(2,2 + dy) ≈ (1 - dy) * A2[2,2] + dy * A2[2,3]
            end
        end
    end

    # AbstractArrays
    makesharedarray(::Type{T}, dims; kwargs...) where {T} = SharedArray{T}(dims; kwargs...)
    function copyshared(A)
        B = makesharedarray(eltype(A), size(A))
        copyto!(B, A)
    end

    for (constructor, copier) in ((interpolate, x->x), (interpolate!, copyshared))
        A2 = makesharedarray(Float64, (N,N), init=A->rand!(A))
        for i = 1:length(A2)
            A2[i] *= 100
        end
        for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
            itp_a = @inferred(constructor(copier(A2), (BSpline(Linear()), BSpline(Quadratic(BC(GT()))))))
            itp_b = @inferred(constructor(copier(A2), (BSpline(Quadratic(BC(GT()))), BSpline(Linear()))))
            if constructor == interpolate!
                @test isa(itp_a.coefs, SharedArray)
                @test isa(itp_b.coefs, SharedArray)
            end
            isfullsize = constructor == interpolate || BC==Periodic
            if isfullsize
                @test @inferred(size(itp_a)) == size(A2)
                @test @inferred(size(itp_b)) == size(A2)
                @test @inferred(axes(itp_a)) == axes(A2)
                @test @inferred(axes(itp_b)) == axes(A2)
            else
                @test @inferred(size(itp_a)) == (N, N-2)
                @test @inferred(size(itp_b)) == (N-2, N)
                @test @inferred(axes(itp_a)) == (1:N, 2:N-1)
                @test @inferred(axes(itp_b)) == (2:N-1, 1:N)
            end

            for j = 2:N-1, i = 2:N-1
                @test ≈(itp_a[i,j],A2[i,j],atol=sqrt(eps(A2[i,j])))
                @test ≈(itp_b[i,j],A2[i,j],atol=sqrt(eps(A2[i,j])))
            end

            for i = 1:10
                dx, dy = rand(), rand()
                @test itp_a(2 + dx,2) ≈ (1 - dx) * A2[2,2] + dx * A2[3,2]
                @test itp_b(2,2 + dy) ≈ (1 - dy) * A2[2,2] + dy * A2[2,3]
            end
        end
    end
end
