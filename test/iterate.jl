using Test
using Interpolations
using Interpolations: KnotRange

# Unit Tests for Base.IteratorEltype and Base.IteratorSize methods (Type Info Only)
@testset "iterate - interface" begin
    import Interpolations.KnotIterator
    # Always have an eltype as we explicitly track via T
    @test Base.IteratorEltype(KnotIterator) == Base.HasEltype()
    @test Base.IteratorEltype(KnotIterator{Int}) == Base.HasEltype()

    # If missing ET type parameter -> SizeUnknown, as could be HasLength or
    # IsInfinite
    @test Base.IteratorSize(KnotIterator) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotIterator{Int}) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotIterator{Int,Flat}) == Base.HasLength()
end

# Iterator Interface for KnotRange
@testset "KnotRange - iterate interface" begin
    import Interpolations.KnotRange
    # Always have an eltype as we explicitly track via T
    @test Base.IteratorEltype(KnotRange) == Base.HasEltype()
    @test Base.IteratorEltype(KnotRange{Int}) == Base.HasEltype()

    # IteratorSize Stored Directly in KnotRange Type
    @test Base.IteratorSize(KnotRange) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotRange{Int}) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotRange{Int,Base.IsInfinite}) == Base.IsInfinite()
    @test Base.IteratorSize(KnotRange{Int,Base.HasLength}) == Base.HasLength()
end

# eltype units tests of KnotIterator / KnotRange
@testset "iterate - eltype" for T ∈ [ Int, Float64, Any ]
    x = convert(Vector{T}, collect(1:5))
    itp = LinearInterpolation(x, x.^2)
    kiter = knots(itp)
    @test Base.IteratorEltype(typeof(kiter)) == Base.HasEltype()
    @test typeof(kiter) <: Interpolations.KnotIterator{T}
    @test eltype(kiter) == T

    krange = knotsbetween(itp; start=1.2, stop=3.4)
    @test Base.IteratorEltype(krange) == Base.HasEltype()
    @test typeof(krange) <: Interpolations.KnotRange{T}
    @test eltype(krange) == T
end

# size/length units tests for KnotIterator
@testset "iterate - size - directional"  begin
    ExtrapBC = [Throw(), Flat(), Periodic(), Reflect()]
    @testset "RevBC: $RevBC, FwdBC: $FwdBC" for FwdBC ∈ ExtrapBC, RevBC ∈ ExtrapBC
        etp = LinearInterpolation(1:10, rand(10),
            extrapolation_bc=((RevBC, FwdBC),)
        )
        kiter = knots(etp)
        @test typeof(kiter) <: Interpolations.KnotIterator
        ktype = KnotIterator{Int, Tuple{typeof(RevBC), typeof(FwdBC)}}
        @test typeof(kiter) == ktype

        if typeof(FwdBC) <: Union{Periodic, Reflect}
            @test Base.IteratorSize(kiter) == Base.IsInfinite()
            @test_throws MethodError length(kiter)
            @test_throws MethodError size(kiter)
        else
            @test Base.IteratorSize(kiter) == Base.HasLength()
            @test length(kiter) == 10
            @test size(kiter) == (10,)
        end
    end
end

# Unit Tests for 1D Interpolations
@testset "iterate - 1d interpolation" begin
    itp = interpolate(1:5, BSpline(Linear()))
    kiter = knots(itp)

    # Single axis -> KnotIterator
    @test typeof(kiter) <: Interpolations.KnotIterator

    # Check size, length and eltype
    @test Base.IteratorSize(kiter) === Base.HasLength()
    @test length(kiter) == 5
    @test size(kiter) == (5,)
    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) <: Int

    # Check contents
    @test collect(kiter) == [1, 2, 3, 4, 5]

    # Sample for loop
    kout = []
    for (k,) ∈ kiter
        push!(kout, k)
    end
    @test kout == collect(kiter)
end

# Unit Test for 2D Interpolations
@testset "iterate - 2d interpolation" begin
    itp = interpolate(rand(3,3), BSpline(Constant()))
    kiter = knots(itp)

    # 2D iterator -> ProductIterator wrapping KnotIterators
    @test typeof(kiter) <: Iterators.ProductIterator
    @test eltype(kiter.iterators) <: Interpolations.KnotIterator

    # Check size, length and eltype
    @test Base.IteratorSize(kiter) == Base.HasShape{2}()
    @test length(kiter) == 9
    @test size(kiter) == (3, 3)
    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) <: Tuple{Int, Int}

    # Collect knots and check again
    k = collect(kiter)
    @test length(k) == 9
    @test size(k) == (3,3)
    @test eltype(k) <: Tuple{Int, Int}
    @test k == Iterators.product(1:3, 1:3) |> collect

    # Sample for loop
    kout = []
    for (kx, ky) ∈ kiter
        push!(kout, (kx, ky))
    end
    @test kout == Iterators.product(1:3, 1:3) |> collect |> vec

end

# Unit Tests for 1D extrapolations
@testset "1D - iterate - $itp" for itp ∈ [ BSpline(Constant()) ]
    x = 1.0:5.0
    itp = interpolate(x.^2, itp)
    @test typeof(itp) <: AbstractInterpolation
    kiter = knots(itp)

    # Checkout construction before proceeding
    @test typeof(kiter) <: Interpolations.KnotIterator
    @test eltype(kiter) <: Int
    @test length(kiter) == 5
    @test size(kiter) == (5,)
    @test collect(kiter) == x

    # Non-repeating Knots -> Should iterate over once and be done
    @testset "extrapolation - $etp" for etp ∈ [ Throw(), Flat(), Line()]
        extrp = extrapolate(itp, etp)
        kiter = knots(extrp)
        @test typeof(kiter) <: Interpolations.KnotIterator
        @test eltype(kiter) <: Int
        @test length(kiter) == 5
        @test size(kiter) == (5,)
        @test collect(kiter) == x
    end

    # Repeating Knots -> Should Iterate indefinitely
    @testset "extrapolation - Periodic" begin
        extrp = extrapolate(itp, Periodic())
        k = knots(extrp)
        @test typeof(k) <: Interpolations.KnotIterator
        @test_throws MethodError length(k)
        @test_throws MethodError size(k)
        k10 = Iterators.take(k, 10) |> collect
        @test k10 ≈ collect(1:10)
        @test extrp.(k10) ≈ vcat(x[1:end-1].^2, x[1:end-1].^2, x[1:2].^2)
    end
    @testset "extrapolation - Reflect" begin
        extrp = extrapolate(itp, Reflect())
        k = knots(extrp)
        @test typeof(k) <: Interpolations.KnotIterator

        # Check length and size throw Method Errors (As Undefined)
        @test Base.IteratorSize(k) == Base.IsInfinite()
        @test_throws MethodError length(k)
        @test_throws MethodError size(k)

        k10 = Iterators.take(k, 10) |> collect
        @test k10 ≈ collect(1:10)
        @test extrp.(k10) ≈ vcat(x.^2, reverse(x[1:end - 1].^2), x[2].^2)
    end
end

# Unit tests for iteration on an uneven grid - Periodic
# Knots should repeat indefinitely with the first and last knots being co-located
@testset "iterate - uneven - Periodic" begin
    x = [1.0, 1.3, 2.4, 3.2, 4.0]
    etp = LinearInterpolation(x, x.^2, extrapolation_bc=Periodic())
    kiter = knots(etp)

    @test typeof(kiter) <: Interpolations.KnotIterator
    @test Base.IteratorSize(kiter) == Base.IsInfinite()
    @test_throws MethodError length(kiter)
    @test_throws MethodError size(kiter)

    k = Iterators.take(kiter, 10) |> collect
    @test k == [1.0, 1.3, 2.4, 3.2, 4.0, 4.3, 5.4, 6.2, 7.0, 7.3]
    @test etp.(k) ≈ vcat(x[1:end-1], x[1:end-1], x[1:2]).^2
end

# Unit tests for iteration on an uneven grid - Reflect
@testset "iterate - uneven - Reflect" begin
    x = [1.0, 1.3, 2.4, 3.2, 4.0]
    etp = LinearInterpolation(x, x.^2, extrapolation_bc=Reflect())
    kiter = knots(etp)

    @test typeof(kiter) <: Interpolations.KnotIterator
    @test Base.IteratorSize(kiter) == Base.IsInfinite()
    @test_throws MethodError length(kiter)
    @test_throws MethodError size(kiter)

    k = Iterators.take(kiter, 10) |> collect
    @test k == [1.0, 1.3, 2.4, 3.2, 4.0, 4.8, 5.6, 6.7, 7.0, 7.3]
    @test etp.(k) ≈ vcat(x, reverse(x[1:end - 1]), x[2]).^2
end

# Test Cases for iteration over the knots of 2D interpolant
#   itersym Iterator under tests (ie. iter = knots(itp))
#   type    The iterator type wrapped by Iterators.ProductIterator
#   expect  One or more iterables for the knots along each dimension
macro test_knots(itersym, type, expect...)

    # 1D vs ND Checks
    if length(expect) > 1
        refExpr = :( Iterators.product($(expect...)) )
        knotType = eval.(expect) .|> eltype |> x -> Tuple{x...}
        typecheckExpr = quote
            @test typeof($itersym) <: Iterators.ProductIterator
            @test typeof($itersym.iterators) <: Tuple{$type, $type}
        end
    else
        refExpr = :( Iterators.flatten($(expect[1])) )
        knotType = eval(expect[1]) |> eltype
        typecheckExpr = :( @test typeof($itersym) <: $type{$knotType} )
    end

    # Assemble component tests into complete testset
    quote
        local $itersym = $(esc(itersym))
        ref = $refExpr
        @testset "typechecks" begin $typecheckExpr end
        @testset "eltype checks" begin
            @test Base.IteratorEltype($itersym) == Base.HasEltype()
            @test eltype($itersym) == $knotType
        end
        @testset "IteratorSize, length and size checks" begin
            @test Base.IteratorSize($itersym) == Base.IteratorSize(ref)
            if Base.IteratorSize(ref) != Base.IsInfinite()
                @test length($itersym) == length(ref)
            end
        end
        @testset "check contents" begin
            if Base.IteratorSize(ref) != Base.IsInfinite()
                kiter = collect($itersym)
                kref = collect(ref)
            else
                kiter = collect(Iterators.take($itersym, 10))
                kref = collect(Iterators.tail($itersym, 10))
            end
            # Test Iterator Results
            @test typeof(kiter) <: AbstractArray
            @test eltype(kiter) == $knotType
            @test length(kiter) == length(kref)
            @test size(kiter) == size(kref)
            @test kiter == kref
        end
    end
end

# Unit tests for 2D iteration with directional boundary conditions that are
# bounded (ie. knots do not repeat indefinitely)
@testset "2D - iteration - bounded - $bc" for bc ∈ [Line(), (Throw(), Line())]
    etp = LinearInterpolation(([1, 2, 3], [1, 2, 3]), rand(3, 3);
        extrapolation_bc=(Throw(), bc)
    )
    kiter = knots(etp)
    @test typeof(kiter) <: Iterators.ProductIterator
    @test eltype(kiter.iterators) <: Interpolations.KnotIterator

    @test Base.IteratorSize(kiter) == Base.HasShape{2}()
    @test length(kiter) == 9
    @test size(kiter) == (3, 3)

    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) == Tuple{Int,Int}

    k = Iterators.take(kiter, 5) |> collect
    @test length(k) == 5
    @test typeof(k) <: Vector{Tuple{Int,Int}}
    kexp = Iterators.product(1:3, 1:3) |> x -> Iterators.take(x, 5) |> collect
    @test k == kexp
end

# Unit tests for 2D iteration with directional boundary conditions that are
# unbounded (ie. knots do repeat indefinitely)
@testset "2D - iteration - Unbounded - $bc" for bc ∈ [Periodic(), (Throw(), Periodic())]
    etp = LinearInterpolation(([1, 2, 3], [1, 2, 3]), rand(3, 3);
        extrapolation_bc=(Line(), bc)
    )
    kiter = knots(etp)
    @test typeof(kiter) <: Iterators.ProductIterator
    @test eltype(kiter.iterators) <: Interpolations.KnotIterator

    # Check length and size throw ArgumentErrors (via Iterators.product)
    @test Base.IteratorSize(kiter) == Base.IsInfinite()
    @test_throws ArgumentError length(kiter)
    @test_throws ArgumentError size(kiter)

    @test Base.IteratorEltype(kiter) == Base.HasEltype()
    @test eltype(kiter) == Tuple{Int,Int}

    k = Iterators.take(kiter, 20) |> collect
    @test length(k) == 20
    @test typeof(k) <: Vector{Tuple{Int,Int}}
    kexp = Iterators.product(1:3, Iterators.countfrom(1)) |>
        x -> Iterators.take(x, 20) |> collect
    @test k == kexp
end

@testset "knotsbetween - interpolate - 1D" begin
    itp = interpolate(rand(10), BSpline(Linear()))

    @testset "start and stop" begin
        krange = knotsbetween(itp; start=2.1, stop=6.2)
        @test_knots krange KnotRange 3:6
    end
    @testset "start" begin
        krange = knotsbetween(itp; start=3.2)
        @test_knots krange KnotRange 4:10
    end
    @testset "stop" begin
        krange = knotsbetween(itp; stop=8.2)
        @test_knots krange KnotRange 1:8
    end
end

@testset "knotsbetween - interpolate - 2D" begin
    itp = interpolate(rand(10, 10), BSpline(Constant()))
    @testset "start and stop" begin
        krange = knotsbetween(itp; start=(2.1, -1.2), stop=(7.2, 7.8))
        @test_knots krange Interpolations.KnotRange 3:7 1:7
    end
    @testset "start" begin
        krange = knotsbetween(itp; start=(-2, 3.2))
        k1, next = iterate(krange)
        @test k1 == (1, 4)
        @test_knots krange Interpolations.KnotRange 1:10 4:10
    end
    @testset "stop" begin
        krange = knotsbetween(itp; stop=(10, 8.2))
        k1, next = iterate(krange)
        @test k1 == (1, 1)
        @test_knots krange Interpolations.KnotRange 1:9 1:8
    end
end

