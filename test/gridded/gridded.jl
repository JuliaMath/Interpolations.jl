using Interpolations, Test
using OffsetArrays: OffsetArray

@testset "Constant Previous at grid knots" begin
    two_knots = interpolate(([1, 2],), [10, 20], Gridded(Constant{Previous}()))
    @test two_knots(1) == 10
    @test two_knots(1.5) == 10
    @test two_knots(2) == 20
    for knots in ([0.0, 1.0, 2.0, 3.0], [-2.0, -0.5, 0.25, 4.0])
        values = [10.0, 20.0, 30.0, 40.0]
        for constructor in (interpolate, interpolate!)
            itp = constructor((knots,), copy(values), Gridded(Constant{Previous}()))
            for i in eachindex(knots)
                @test itp(knots[i]) == values[i]
                if i > firstindex(knots)
                    @test itp(prevfloat(knots[i])) == values[i-1]
                end
                if i < lastindex(knots)
                    @test itp(nextfloat(knots[i])) == values[i]
                end
            end
            @test @inferred(itp(last(knots))) == last(values)
            midpoints = (knots[1:end-1] .+ knots[2:end]) ./ 2
            for i in eachindex(midpoints)
                @test itp(midpoints[i]) == values[i]
            end
            xs = [knots[1], midpoints[1], knots[2], midpoints[2],
                  knots[3], midpoints[3], knots[4]]
            expected = [10.0, 10.0, 20.0, 20.0, 30.0, 30.0, 40.0]
            @test itp(xs) == expected
            @test itp(reverse(xs)) == reverse(expected)
            @test_throws BoundsError itp(prevfloat(first(knots)))
            @test_throws BoundsError itp(nextfloat(last(knots)))
            flat = extrapolate(itp, Flat())
            @test flat(first(knots) - 1) == first(values)
            @test flat(last(knots) + 1) == last(values)
        end

        next = interpolate((knots,), values, Gridded(Constant{Next}()))
        nearest = interpolate((knots,), values, Gridded(Constant()))
        linear = interpolate((knots,), values, Gridded(Linear()))
        for i in 1:length(knots)-1
            midpoint = (knots[i] + knots[i+1]) / 2
            @test next(midpoint) == values[i+1]
            @test nearest(midpoint) == values[i+1] # ties go to the next knot
            @test linear(midpoint) == (values[i] + values[i+1]) / 2
        end
        @test nearest(knots) == values
        @test linear(knots) == values
    end

    xs, ys = [1.0, 2.0, 4.0], [0.0, 3.0]
    values = [11.0 12.0; 21.0 22.0; 31.0 32.0]
    itp = interpolate((xs, ys), values, Gridded(Constant{Previous}()))
    for j in eachindex(ys), i in eachindex(xs)
        @test itp(xs[i], ys[j]) == values[i,j]
    end
    @test itp([1.0, 1.5, 2.0, 3.0, 4.0], [0.0, 1.5, 3.0]) ==
        [11.0 11.0 12.0; 11.0 11.0 12.0; 21.0 21.0 22.0;
         21.0 21.0 22.0; 31.0 31.0 32.0]
    mixed = interpolate((xs, ys), values,
        (Gridded(Constant{Previous}()), Gridded(Linear())))
    @test mixed(2.0, 1.5) == 21.5
    @test mixed(4.0, 1.5) == 31.5

    mixed_knots = interpolate((2:4, [0, 1]), values, Gridded(Constant{Previous}()))
    @test mixed_knots(2, 0) == values[1, 1]
    @test mixed_knots(3, 1) == values[2, 2]
    @test mixed_knots(4, 1) == values[3, 2]
    @test mixed_knots(2.5, 1) == values[1, 2]
    @test mixed_knots(2:4, [0, 1]) == values

    offset_knots = OffsetArray([1.0, 2.0, 4.0], -1:1)
    offset_values = OffsetArray([10.0, 20.0, 40.0], -1:1)
    offset_itp = interpolate((offset_knots,), offset_values, Gridded(Constant{Previous}()))
    for i in eachindex(offset_knots)
        @test offset_itp(offset_knots[i]) == offset_values[i]
    end
    @test offset_itp(3.0) == 20.0
end

@testset "LinearTests" begin
    front(r::AbstractUnitRange) = first(r):last(r)-1
    front(r::AbstractRange) = range(first(r), step=step(r), length=length(r)-1)

    for D in (Constant, Linear)
        for constructor in (interpolate, interpolate!)
            isinplace = constructor == interpolate!
            ## 1D
            a = rand(5)
            knots = (range(1, stop=length(a), length=length(a)),)
            itp = @inferred(constructor(knots, a, Gridded(D())))
            @inferred(itp(2))
            @inferred(itp(CartesianIndex(2)))
            for i = 1:length(a)
                @test itp(i) ≈ a[i]
                @test itp(CartesianIndex(i)) ≈ a[i]
            end
            @inferred(itp(knots...))
            @test itp(knots...) ≈ a
            # compare scalar indexing and vector indexing
            x = front(knots[1] .+ 0.1)
            v = itp(x)
            for i = 1:length(x)
                @test v[i] ≈ itp(x[i])
            end
            x = [2.3,2.2]   # non-increasing order
            v = itp(x)
            for i = 1:length(x)
                @test v[i] ≈ itp(x[i])
            end
            # compare against BSpline
            itpb = @inferred(constructor(a, BSpline(D())))
            for x in range(1.1, stop=4.9, length=101)
                @test itp(x) ≈ itpb(x)
            end

            knots = (range(0.0, stop=1.0, length=length(a)),)
            itp = @inferred(constructor(knots, a, Gridded(D())))
            @test itp([0.1, 0.2, 0.3]) == [itp(0.1), itp(0.2), itp(0.3)]

            ## 2D
            A = rand(6,5)
            knots = (range(1, stop=size(A,1), length=size(A,1)), range(1, stop=size(A,2), length=size(A,2)))
            itp = @inferred(constructor(knots, A, Gridded(D())))
            isinplace && @test parent(itp) === A
            @inferred(itp(2, 2))
            @inferred(itp(CartesianIndex((2,2))))
            for j = 2:size(A,2)-1, i = 2:size(A,1)-1
                @test itp(i,j) ≈ A[i,j]
                @test itp(CartesianIndex((i,j))) ≈ A[i,j]
            end
            @test itp(knots...) ≈ A
            @inferred(itp(knots...))
            # compare scalar indexing and vector indexing
            x, y = front(knots[1] .+ 0.1), front(knots[2] .+ 0.6)
            v = itp(x,y)
            for j = 1:length(y), i = 1:length(x)
                @test v[i,j] ≈ itp(x[i],y[j])
            end
            # check the fallback vector indexing
            x = [2.3,2.2]   # non-increasing order
            y = [3.5,2.8]
            v = itp(x,y)
            for j = 1:length(y), i = 1:length(x)
                @test v[i,j] ≈ itp(x[i],y[j])
            end
            # compare against BSpline
            itpb = @inferred(constructor(A, BSpline(D())))
            for x in range(1.1, stop=5.9, length=101), y in range(1.1, stop=4.9, length=101)
                @test itp(x,y) ≈ itpb(x,y)
            end

            A = rand(8,20)
            knots = ([x^2 for x = 1:8], [0.2y for y = 1:20])
            itp = constructor(knots, A, Gridded(D()))
            @test itp(4,1.2) ≈ A[2,6]
        end
    end

    # issue #248
    itp = interpolate((1001:1005,), 1:5, Gridded(Linear()))
    @test itp(1002) ≈ 2

    # issue #255
    itp = interpolate(([1,2,3],), [1.0f0, 2.0f0, 1.0f0], Gridded(Linear()))
    @test itp(2) === 2.0f0
    @test itp(2.0) === 2.0
    @test itp(2.0f0) === 2.0f0

    # trailing 1s, issue #301
    @test itp(1.8, 1) == itp(1.8)

    # issue #395
    knots = ([0.0,0.0,1.0],) # not unique
    # Further discussion determined this should be a warning rather than an error
    # https://github.com/JuliaMath/Interpolations.jl/commit/318ebc88ca1fc084754f3c741266537f901a3310
    knots_not_unique_warning = "Duplicated knots were deduplicated. Use Interpolations.deduplicate_knots!(knots) explicitly to avoid this warning."
    @test_logs (:warn, knots_not_unique_warning) interpolate(knots, [1.0, 1.1, 2.0], Gridded(Linear()))

    # knot deduplication, issue #467, PR #468
    duplicated_knots = [1.0, 1.0, 1.0, nextfloat(1.0), nextfloat(1.0)]
    successive_knots_warning = "Successive repeated knots detected. Consider using `move_knots` keyword to Interpolations.deduplicate_knots!"
    @test_logs (:warn, successive_knots_warning) Interpolations.deduplicate_knots!(duplicated_knots)
    @test allunique(duplicated_knots)

    # inplace gridded interpolation, issue #495
    knots = ([0., 1.],)
    y = view([1., 2.], :)
    f1 = interpolate(knots, y, Gridded(Linear()))
    f2 = interpolate!(knots, y, Gridded(Linear()))
    y .= 0
    @test f1(0.5) ≈ 1.5
    @test f2(0.5) ≈ 0
end
