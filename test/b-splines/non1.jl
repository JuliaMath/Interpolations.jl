module Non1Tests

using Interpolations, OffsetArrays, AxisAlgorithms, Test

# At present, for a particular type of non-1 array you need to specialize this function
function AxisAlgorithms.A_ldiv_B_md!(dest::OffsetArray, F, src::OffsetArray, dim::Integer, b::AbstractVector)
    indsdim = axes(parent(src), dim)
    indsF = axes(F)[2]
    if indsF == indsdim
        return A_ldiv_B_md!(parent(dest), F, parent(src), dim, b)
    end
    throw(DimensionMismatch("indices $(axes(parent(src))) do not match $(axes(F))"))
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
        @test @inferred(axes(itp1)) === axes(A1)

        # test that we reproduce the values at on-grid points
        for x = inds
            @test itp1[x] == itp1(x) ≈ f1(x)
        end

        itp2 = @inferred(constructor(copier(A2), BSpline(O()), GT()))
        @test @inferred(axes(itp2)) === axes(A2)
        for j = yinds, i = xinds
            @test itp2[i,j] == itp2(i,j) ≈ A2[i,j]
        end
    end

    for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
        itp1 = @inferred(constructor(copier(A1), BSpline(Quadratic(BC())), GT()))
        isfullsize = constructor == interpolate || BC==Periodic
        if isfullsize
            @test @inferred(size(itp1)) == size(A1)
            @test @inferred(axes(itp1)) == axes(A1)
        else
            @test @inferred(size(itp1)) == (length(Base.Slice(first(inds)+1:last(inds)-1)),)
            @test @inferred(axes(itp1)) == (Base.Slice(first(inds)+1:last(inds)-1),)
        end

        for x in eachindex(itp1)
            @test itp1[x] == itp1(x) ≈ f1(x)
        end

        itp2 = @inferred(constructor(copier(A2), BSpline(Quadratic(BC())), GT()))
        if isfullsize
            @test @inferred(size(itp2)) == size(A2)
            @test @inferred(axes(itp2)) == axes(A2)
        else
            @test @inferred(size(itp2)) == (29, 8)
            @test @inferred(axes(itp2)) == (Base.Slice(-1:27), Base.Slice(1:8))
        end
        for x in eachindex(itp2)
            @test itp2[x] == itp2(x) ≈ A2[x]
        end
    end

    for BC in (Flat,Line,Free,Periodic), GT in (OnGrid, OnCell)
        itp1 = @inferred(constructor(copier(A1), BSpline(Cubic(BC())), GT()))
        isfullsize = constructor == interpolate || BC==Periodic
        if isfullsize
            @test @inferred(size(itp1)) == size(A1)
            @test @inferred(axes(itp1)) == axes(A1)
        else
            @test @inferred(size(itp1)) == (length(Base.Slice(first(inds)+1:last(inds)-1)),)
            @test @inferred(axes(itp1)) == (Base.Slice(first(inds)+1:last(inds)-1),)
        end

        # test that we reproduce the values at on-grid points
        for x = eachindex(itp1)
            @test itp1[x] == itp1(x) ≈ f1(x)
        end

        itp2 = @inferred(constructor(copier(A2), BSpline(Cubic(BC())), GT()))
        if isfullsize
            @test @inferred(size(itp2)) == size(A2)
            @test @inferred(axes(itp2)) == axes(A2)
        else
            @test @inferred(size(itp2)) == (29, 8)
            @test @inferred(axes(itp2)) == (Base.Slice(-1:27), Base.Slice(1:8))
        end
        for x in eachindex(itp2)
            @test_skip itp2[x] == itp2(x) ≈ A2[x]
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
        @test itp2[i,j] == itp2(i,j) ≈ A2[i,j]
    end
end

end
