module Non1Tests

using Interpolations, OffsetArrays, AxisAlgorithms, Base.Test

# At present, for a particular type of non-1 array you need to specialize this function
function AxisAlgorithms.A_ldiv_B_md!(dest::OffsetArray, F, src::OffsetArray, dim::Integer, b::AbstractVector)
    indsdim = indices(parent(src), dim)
    indsF = indices(F)[2]
    if indsF == indsdim
        return A_ldiv_B_md!(parent(dest), F, parent(src), dim, b)
    end
    throw(DimensionMismatch("indices $(indices(parent(src))) do not match $(indices(F))"))
end

for (constructor, copier) in ((interpolate, x->x), (interpolate!, copy))
    f1(x) = sin((x-3)*2pi/9 - 1)
    inds = -3:6
    A1 = OffsetArray(Float64[f1(x) for x in inds], inds)

    f2(x,y) = sin(x/10)*cos(y/6) + 0.1
    xinds, yinds = -2:28,0:9
    A2 = OffsetArray(Float64[f2(x,y) for x in xinds, y in yinds], xinds, yinds)

    for GT in (OnGrid, OnCell), O in (Constant, Linear)
        itp1 = @inferred(constructor(copier(A1), BSpline(O()), GT()))
        @test @inferred(indices(itp1)) === indices(A1)

        # test that we reproduce the values at on-grid points
        for x = inds
            @test itp1[x] ≈ f1(x)
        end

        itp2 = @inferred(constructor(copier(A2), BSpline(O()), GT()))
        @test @inferred(indices(itp2)) === indices(A2)
        for j = yinds, i = xinds
            @test itp2[i,j] ≈ A2[i,j]
        end
    end

    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        itp1 = @inferred(constructor(copier(A1), BSpline(Quadratic(BC())), GT()))
        @test @inferred(indices(itp1)) === indices(A1)

        # test that we reproduce the values at on-grid points
        inset = constructor == interpolate!
        for x = first(inds)+inset:last(inds)-inset
            @test itp1[x] ≈ f1(x)
        end

        itp2 = @inferred(constructor(copier(A2), BSpline(Quadratic(BC())), GT()))
        @test @inferred(indices(itp2)) === indices(A2)
        for j = first(yinds)+inset:last(yinds)-inset, i = first(xinds)+inset:last(xinds)-inset
            @test itp2[i,j] ≈ A2[i,j]
        end
    end

    for BC in (Flat,Line,Free,Periodic), GT in (OnGrid, OnCell)
        itp1 = @inferred(constructor(copier(A1), BSpline(Cubic(BC())), GT()))
        @test @inferred(indices(itp1)) === indices(A1)

        # test that we reproduce the values at on-grid points
        inset = constructor == interpolate!
        for x = first(inds)+inset:last(inds)-inset
            @test itp1[x] ≈ f1(x)
        end

        itp2 = @inferred(constructor(copier(A2), BSpline(Cubic(BC())), GT()))
        @test @inferred(indices(itp2)) === indices(A2)
        for j = first(yinds)+inset:last(yinds)-inset, i = first(xinds)+inset:last(xinds)-inset
            @test itp2[i,j] ≈ A2[i,j]
        end
    end
end

let
    f(x) = sin((x-3)*2pi/9 - 1)
    inds = -7:2
    A = OffsetArray(Float64[f(x) for x in inds], inds)
    itp1 = interpolate!(copy(A), BSpline(Quadratic(InPlace())), OnCell())
    for i in inds
        @test itp1[i] ≈ A[i]
    end

    f(x,y) = sin(x/10)*cos(y/6) + 0.1
    xinds, yinds = -2:28,0:9
    A2 = OffsetArray(Float64[f(x,y) for x in xinds, y in yinds], xinds, yinds)
    itp2 = interpolate!(copy(A2), BSpline(Quadratic(InPlace())), OnCell())
    for j = yinds, i = xinds
        @test itp2[i,j] ≈ A2[i,j]
    end
end

end
