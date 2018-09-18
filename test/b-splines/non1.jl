using AxisAlgorithms, OffsetArrays

@testset "Unconventional axes" begin
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
        isinplace = constructor == interpolate!
        f1(x) = sin((x-3)*2pi/9 - 1)
        inds = -3:6
        A1 = OffsetArray(Float64[f1(x) for x in inds], inds)

        f2(x,y) = sin(x/10)*cos(y/6) + 0.1
        xinds, yinds = -2:28,0:9
        A2 = OffsetArray(Float64[f2(x,y) for x in xinds, y in yinds], xinds, yinds)

        for O in (Constant, Linear)
            itp1 = @inferred(constructor(copier(A1), BSpline(O())))
            check_axes(itp1, A1, isinplace)
            check_inbounds_values(itp1, A1)
            check_oob(itp1)
            can_eval_near_boundaries(itp1)

            itp2 = @inferred(constructor(copier(A2), BSpline(O())))
            check_axes(itp2, A2, isinplace)
            check_inbounds_values(itp2, A2)
            check_oob(itp2)
            can_eval_near_boundaries(itp2)
        end

        for BC in (Flat,Line,Free,Periodic,Reflect,Natural), GT in (OnGrid, OnCell)
            itp1 = @inferred(constructor(copier(A1), BSpline(Quadratic(BC(GT())))))
            check_axes(itp1, A1, isinplace)
            check_inbounds_values(itp1, A1)
            check_oob(itp1)
            can_eval_near_boundaries(itp1)

            itp2 = @inferred(constructor(copier(A2), BSpline(Quadratic(BC(GT())))))
            check_axes(itp2, A2, isinplace)
            check_inbounds_values(itp2, A2)
            check_oob(itp2)
            can_eval_near_boundaries(itp2)
        end

        for BC in (Flat,Line,Free,Periodic), GT in (OnGrid, OnCell)
            itp1 = @inferred(constructor(copier(A1), BSpline(Cubic(BC(GT())))))
            check_axes(itp1, A1, isinplace)
            check_inbounds_values(itp1, A1)
            check_oob(itp1)
            can_eval_near_boundaries(itp1)

            itp2 = @inferred(constructor(copier(A2), BSpline(Cubic(BC(GT())))))
            check_axes(itp2, A2, isinplace)
            check_inbounds_values(itp2, A2)
            check_oob(itp2)
            can_eval_near_boundaries(itp2)
        end
    end

    let
        f(x) = sin((x-3)*2pi/9 - 1)
        inds = -7:2
        A = OffsetArray(Float64[f(x) for x in inds], inds)
        itp1 = interpolate!(copy(A), BSpline(Quadratic(InPlace(OnCell()))))
        check_axes(itp1, A)
        check_inbounds_values(itp1, A)
        check_oob(itp1)
        can_eval_near_boundaries(itp1)

        f(x,y) = sin(x/10)*cos(y/6) + 0.1
        xinds, yinds = -2:28,0:9
        A2 = OffsetArray(Float64[f(x,y) for x in xinds, y in yinds], xinds, yinds)
        itp2 = interpolate!(copy(A2), BSpline(Quadratic(InPlace(OnCell()))))
        check_axes(itp2, A2)
        check_inbounds_values(itp2, A2)
        check_oob(itp2)
        can_eval_near_boundaries(itp2)
    end
end
