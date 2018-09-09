using Interpolations, Test, ForwardDiff

@testset "Issues" begin
    @testset "issue 34" begin
        A = rand(1:20, 100, 100)

        # In #34, this incantation throws
        itp = interpolate(A, BSpline(Quadratic(Flat(OnCell()))))
        # Sanity check that not only don't throw, but actually interpolate
        for i in 1:size(A,1), j in 1:size(A,2)
            @test itp[i,j] ≈ A[i,j]
        end
    end

    @testset "issue 129" begin
        xy = [3z+2y for z in range(0.0,stop=1,length=10),y in 0:0.2:1]
        itp = interpolate((1:10,1:6),xy,(Gridded(Linear()),NoInterp()))
        @test itp(3.3,1:6) == [itp(3.3,i) for i = 1:6]
        sitp = scale(itp ,range(0.0,stop=1.0,length=10),1:6)
        @test sitp(0.8,1:6) == [sitp(0.8,i) for i = 1:6]
    end

    @testset "issue 151" begin
        V = zeros(10,10,10,10)
        interpV = interpolate(V, BSpline(Cubic(Line(OnGrid()))))
        @test ndims(interpV) == 4
    end

    @testset "issue 158" begin
        A_x = 1.:2.:40.
        A = [log(x) for x in A_x]
        itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
        sitp = scale(itp, A_x)
        @test_throws BoundsError sitp(42.0)
        @test_throws BoundsError sitp([3.0,42.0])
    end

    @testset "issue 165" begin
        V0 = Array{Float64, 2}(undef, 5, 6)
        yGrid = range(-1.0,stop= 1.0,length= 5)
        bGrid = range(2.0, stop=3.0, length=6)
        iV = scale(interpolate(V0, BSpline(Cubic(Line(OnCell())))), yGrid, bGrid)
        @test iV isa Interpolations.AbstractInterpolation
    end

    @testset "issue 191" begin
        function eff(p, thing)
            out = p[1] *  p[2]
            return out
        end
        function GradientWrapper(func, p, thing)
            f(x) = func(x, thing)
            return ForwardDiff.gradient(f,p)
        end
        io = IOBuffer()
        println(io, "Start.")
        n = 40
        z = [i * j for i in 1:n, j in 1:n]
        itp = interpolate(z, BSpline(Cubic(Natural(OnGrid()))))
        p = [5.0, 6.0]
        test = [itp, itp]
        println(io, "Using type $(typeof(test))")
        println(io, "Calling the derivative directly:")
        f(x) = eff(x,test)
        println(io, ForwardDiff.gradient(f, p))
        println(io, "Calling the wrapper:")
        println(io, GradientWrapper(eff, p, test))
        test = Array{Any}(undef, 2)
        test[1] = itp
        test[2] = itp
        println(io, "Using type $(typeof(test))")
        println(io, "Calling the derivative directly:")
        g(x) = eff(x, test)
        println(io, ForwardDiff.gradient(g, p))
        println(io, "Calling the wrapper:")
        println(io, GradientWrapper(eff, p, test))
        test = Array{AbstractInterpolation}(undef, 2)
        test[1] = itp
        test[2] = itp
        println(io, "Using type $(typeof(test))")
        println(io, "Calling the derivative directly:")
        h(x) = eff(x, test)
        println(io, ForwardDiff.gradient(h, p))
        println(io, "Calling the wrapper -- this fails:")
        println(io, GradientWrapper(eff, p, test))
        println(io, "Done.")
        str = String(take!(io))
        @show str
        @test endswith(str, "Done.\n")
    end

    @testset "issue 200" begin
        grid = range(0, stop=2.0, length=10)
        vals = [-exp(-x) for x in grid]
        itp = interpolate(vals, BSpline(Cubic(Line(OnGrid()))))
        sitp = scale(itp, grid)
        s = sitp(0.0)
        @test sitp([0.0, 0.0]) == [s, s]
    end

    @testset "issue 202" begin
        xs = 1:10
        itp = scale(interpolate(zeros(10), BSpline(Cubic(Periodic(OnGrid())))), xs)
        extp = extrapolate(itp, Periodic())
        @test Interpolations.gradient(extp, 2.) == [0.0]
        @test_throws ErrorException Interpolations.gradient(extp, 2.0, 2.0)

        # The issue title says 3d, so let's try a 3d case
        A = (0.2.*(1:3)) .* (0.3.*(1:3))' .* reshape(0.4.*(1:3), 1, 1, 3)
        itp = interpolate(A, BSpline(Linear()))
        etp = extrapolate(itp, Line())
        @test Interpolations.gradient(etp, 2, 2, 2) ≈ [0.2*0.6*0.8, 0.4*0.3*0.8, 0.4*0.6*0.4]
        @test Interpolations.gradient(etp, 10, 2, 2) ≈ [0.2*0.6*0.8, 0.6*0.3*0.8, 0.6*0.6*0.4]
        @test Interpolations.gradient(etp, 2, 10, 2) ≈ [0.2*0.9*0.8, 0.4*0.3*0.8, 0.4*0.9*0.4]
        @test Interpolations.gradient(etp, 2, 2, 10) ≈ [0.2*0.6*1.2, 0.4*0.3*1.2, 0.4*0.6*0.4]
    end

    @testset "issue 213" begin
        function build_itp()
            A_x = 1.0:2.0:40.0
            A = A_x.^2
            knots = (A_x,)
            interpolate(knots, A, Gridded(Linear()))
        end
        function itp_test(x::AbstractVector)
            itp = build_itp()
            return itp(x...)
        end

        j = ForwardDiff.gradient(itp_test, [2.0])
        @test j ≈ Interpolations.gradient(build_itp(), 2.0)
    end
end
