@testset "Regularization" begin
    for (constructor, copier) in ((interpolate, identity), (interpolate!, copy))
        isinplace = constructor == interpolate!
        f0(x) = sin((x-3)*2pi/9 - 1)
        f1(x) = 1.0 + 0.1*x + 0.01*x^2 + 0.001*x^3

        xmax = 10
        xi = range(0, stop=xmax, length=11)
        xfine = range(0, stop=xmax, length=101)

        A0 = f0.(xi)
        A1 = f1.(xi)

        for (yi, f) in ((A0, f0), (A1, f1)), λ in (0.01, 100), k in (1, 2, 3, 4)
            it = @inferred(interpolate(yi, BSpline(Cubic(Line(OnCell()))), λ, k))
            itt = @inferred(scale(it, xi))
            
            if f === f0
                if λ == 0.01
                    # test that interpolation of the interpolating points is close to expected data
                    for x in [2, 3, 6, 7]
                        @test f(x) ≈ itt(x) atol=0.01
                    end
                elseif λ == 100
                    # test that interpolation of the interpolating points is far from expected data
                    for x in [2, 3, 6, 7]
                        @test !(isapprox(f(x), itt(x); atol=0.1))
                    end
                end
            
            elseif f === f1
                if k==4
                    for x in xfine
                        # test that interpolation is close to expected data when
                        # the smoothing order is higher than the polynomial order
                        @test f(x) ≈ itt(x) atol=0.1
                    end
                end
            end
        end

        # Check fallback to no regularization for λ=0
        itp1  = @inferred(constructor(copier(A0), BSpline(Cubic(Line(OnGrid()))), 0, 1))
        itp1p = @inferred(constructor(copier(A0), BSpline(Cubic(Line(OnGrid())))))
        @test itp1 == itp1p


        f2(x, y) = sin(x/10)*cos(y/6)
        xmax2, ymax2 = 30, 10
        A2 = Float64[f2(x, y) for x in 1:xmax2, y in 1:ymax2]

        # Check fallback to no regularization for non Vectors
        itp2  = @inferred(constructor(copier(A2), BSpline(Cubic(Line(OnGrid()))), 1, 1))
        itp2p = @inferred(constructor(copier(A2), BSpline(Cubic(Line(OnGrid())))))
        @test itp2 == itp2p
    end
end
