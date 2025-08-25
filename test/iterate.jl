using Test
using Interpolations
using Interpolations: KnotRange

"""
    @test_knots kiter T expect

Generate test sets for testing the knot iterator `kiter`. Following checks are
performed:
- `typeof(kiter)` is `T` if 1D or `ProductIterator` wrapping `T` if ND
- `IteratorEltype` and `eltype` behavior matches expected
- `IteratorSize`, `length` and `size` are matches expected
- Iterated contents match (Checks all if bounded, otherwise the first 100)
"""
macro test_knots(itersym, type, expect...)
    # 1D vs ND Checks
    if length(expect) > 1
        refExpr = :( Iterators.product($(expect...)) )
        knotTypeExpr = :( Tuple{ eltype.(($(expect...),))... } )
        typecheckExpr = quote
            @test typeof($itersym) <: Iterators.ProductIterator
            @test typeof($itersym.iterators) <: Tuple{$type, $type}
        end
        sizecheckExpr = :( @test size($itersym) == size(ref) )

        # Expect HasShape{N} for ND Iterators
        IteratorSizeExpr = :( Base.HasShape{$(length(expect))}() )

        # Content Check for 1D array of tuples
        IterEqCheck = :( all(map( (x,y) -> (all∘map)(≈, x, y), kiter, kref )))
    else
        refExpr = expect[1]
        knotTypeExpr = :( eltype($refExpr) )
        typecheckExpr = :( @test (<:)(typeof($itersym), $type) )

        # Check size using length to avoid issues with reference not defining a
        # suitable size method
        sizecheckExpr = :( @test size($itersym) == (length(ref),) )

        # Expect a IteratorSize of HasLength for 1D iterators
        IteratorSizeExpr = :( Base.HasLength() )

        # Content Check for 1D arrays
        IterEqCheck = :( kiter ≈ kref )
    end

    # Check Type-stability using inferred
    type_stability_checks = quote
        @testset "type stability" begin
            iteratetype = Union{Nothing, knotType}
            y = @inferred iteratetype iterate($itersym)
            if y !== nothing
                @test_nowarn @inferred iteratetype iterate($itersym, y[2])
            end
            isBounded && @test_nowarn @inferred Integer length($itersym)
            @test_nowarn @inferred Type eltype($itersym)
        end
    end

    # Build test_knots tests
    quote
        local $itersym = $(esc(itersym))    # Local knot iterator for testing
        knotType = $(esc(knotTypeExpr))     # Insert Expr to get knot eltype
        ref = $(esc(refExpr))               # Build reference iterator

        # If the reference size is "SizeUnknown" treat as IsInfinite
        # Some iterators are quick to assume SizeUnknown
        isBounded = Base.IteratorSize(ref) != Base.IsInfinite() &&
                    Base.IteratorSize(ref) != Base.SizeUnknown()

        # Insert Type checks defined earlier
        @testset "typechecks" begin $typecheckExpr end

        # Check Eltype matches reference iterator
        # Relax type check (ie. <: instead of ==) as some iterators are quick to
        # assume "Any" -> Prevents errors caused by ref's type being off
        # This is okay, as we check that collected eltypes match exactly later
        @testset "eltype checks" begin
            @test Base.IteratorEltype($itersym) == Base.HasEltype()
            @test eltype($itersym) <: knotType
        end

        # Check that IteratorSize matches reference and iff bounded check that
        # length matches as well
        # Size check is defined by $sizecheckExpr to avoid issues with Undefined
        # size methods for 1D reference iterators
        @testset "IteratorSize, length and size checks" begin
            if isBounded
                @test Base.IteratorSize($itersym) == $IteratorSizeExpr
                @test length($itersym) == length(ref)
                $sizecheckExpr
            else
                @test Base.IteratorSize($itersym) == Base.IsInfinite()
            end
        end

        # Check that the contents of the knot iterator matches the reference
        # iterator
        @testset "check contents" begin
            if isBounded
                kiter = collect($itersym)
                kref = collect(ref)
            else
                kiter = collect(Iterators.take($itersym, 100))
                kref = collect(Iterators.take(ref, 100))
            end
            # Test Iterator Results
            @test typeof(kiter) <: AbstractArray
            @test eltype(kiter) == eltype(kref)
            @test length(kiter) == length(kref)
            @test size(kiter) == size(kref)
            @test $IterEqCheck
        end

        # Type stability -> Check that methods are type suitable
        $type_stability_checks
    end
end

"""
    knot_ref(seq, etp)
    knot_ref(seq, etp; start, stop)

Dumb but correct construction of a reference knot iterator for a boundary condition
of etp.

If `start` is specified, will skip knots until `start < k` using Iterator.filter
If `stop` is specified, will iterate until `stop <= k` and return an array
"""
knot_ref(seq, etp; start=nothing, stop=nothing) = knot_ref(seq, etp, start, stop)
function knot_ref(seq, etp, start::Number, ::Nothing)
    # Drop items less than start -> Dispatch to generate sequence
    Iterators.filter(x -> start < x, knot_ref(seq, etp))
end
function knot_ref(seq, etp, start, stop::Number)
    # Collect items until `k < stop` is false -> Then return ref
    # Dispatch to handle remaining arguments
    ref = Vector{eltype(seq)}()
    for k ∈ knot_ref(seq, etp, start, nothing)
        if k < stop
            push!(ref, k)
        else
            break
        end
    end
    ref
end

knot_ref(seq, etp, ::Nothing, ::Nothing) = knot_ref(seq, etp)

knot_ref(seq, ::Interpolations.BoundaryCondition) = seq

function knot_ref(seq, ::Periodic)
    Δ = diff(seq)
    iter = Iterators.flatten((first(seq), Iterators.cycle(Δ)))
    Iterators.accumulate(+, iter)
end

function knot_ref(seq, ::Reflect)
    Δ = diff(seq)
    Δ = vcat(Δ, reverse(Δ))
    iter = Iterators.flatten((first(seq), Iterators.cycle(Δ)))
    Iterators.accumulate(+, iter)
end

@testset "knot_ref" begin
    x = [1.0, 1.5, 1.75, 2.0]
    # Non-Repeating Knots
    @testset "knot_ref for $etp" for etp ∈ [Line(), Flat(), Throw()]
        @test knot_ref(x, etp) == x
        @test knot_ref(x, etp; start=1.2) |> collect == x[2:end]
        @test knot_ref(x, etp; stop=1.8) |> collect == x[1:end-1]
        @test knot_ref(x, etp; start=1.2, stop=1.8) |> collect == x[2:end-1]
    end
    # Check knot_ref for Periodic against expect knots
    @testset "knot_ref for Periodic" begin
        ref = [
            1.0, 1.5, 1.75, 2.0, 2.5, 2.75, 3.0, 3.5, 3.75, 4.0, 4.5, 4.75,
            5.0, 5.5, 5.75, 6.0, 6.5, 6.75, 7.0, 7.5, 7.75, 8.0, 8.5, 8.75,
            9.0, 9.5, 9.75, 10.0, 10.5, 10.75, 11.0, 11.5, 11.75, 12.0
        ]
        @test knot_ref(x, Periodic(), start=2) |> x -> Iterators.take(x, 10) |> collect == ref[5:14]
        @test knot_ref(x, Periodic(); start=1.0, stop=12.0) == ref[2:end-1]
        @test knot_ref(x, Periodic(); start=2.3, stop=11.0) == ref[2.3 .< ref .< 11.0]
   end
    # Check knot_ref for Reflect against expect knots
    @testset "knot_ref for Reflect" begin
        ref = [
            1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 3.75, 4.0, 4.25, 4.5,
            5.0, 5.5, 5.75, 6.0, 6.25, 6.5, 7.0, 7.5, 7.75, 8.0, 8.25, 8.5,
            9.0, 9.5, 9.75, 10.0, 10.25, 10.5, 11.0, 11.5, 11.75, 12.0
        ]
        @test knot_ref(x, Reflect(), start=2) |> x -> Iterators.take(x, 10) |> collect == ref[5:14]
        @test knot_ref(x, Reflect(); start=1.0, stop=12.0) == ref[2:end-1]
        @test knot_ref(x, Reflect(); start=2.3, stop=11.0) == ref[2.3 .< ref .< 11.0]
    end
end

# Unit Tests for Base.IteratorEltype and Base.IteratorSize methods (Type Info Only)
@testset "iterate - interface" begin
    import Interpolations.KnotIterator
    # Eletype is known iff type parameter is provided
    @test Base.IteratorEltype(KnotIterator) == Base.EltypeUnknown()
    @test eltype(KnotIterator) == Any
    @test Base.IteratorEltype(KnotIterator{Int}) == Base.HasEltype()
    @test eltype(KnotIterator{Int}) == Int
    @test Base.IteratorEltype(KnotIterator{Int,Flat}) == Base.HasEltype()
    @test eltype(KnotIterator{Int}) == Int

    # If missing ET type parameter -> SizeUnknown, as could be HasLength or
    # IsInfinite
    @test Base.IteratorSize(KnotIterator) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotIterator{Int}) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotIterator{Int,Flat}) == Base.HasLength()
end

# Iterator Interface for KnotRange
@testset "KnotRange - iterate interface" begin
    import Interpolations.KnotRange
    # Eletype is known iff type parameter is provided
    @test Base.IteratorEltype(KnotRange) == Base.EltypeUnknown()
    @test eltype(KnotRange) == Any
    @test Base.IteratorEltype(KnotRange{Int}) == Base.HasEltype()
    @test eltype(KnotRange{Int}) == Int
    @test Base.IteratorEltype(KnotRange{Int,Base.UnitRange}) == Base.HasEltype()
    @test eltype(KnotRange{Int,Base.UnitRange}) == Int
    @test Base.IteratorEltype(KnotRange{Int,Iterators.Count}) == Base.HasEltype()
    @test eltype(KnotRange{Int,Iterators.Count}) == Int

    # If missing Range type parameter -> SizeUnknown, as could be HasLength or
    # IsInfinite
    @test Base.IteratorSize(KnotRange) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotRange{Int}) == Base.SizeUnknown()
    @test Base.IteratorSize(KnotRange{Int,Base.UnitRange}) == Base.HasLength()
    @test Base.IteratorSize(KnotRange{Int,Iterators.Count}) == Base.IsInfinite()
end

# eltype units tests of KnotIterator / KnotRange
@testset "iterate - eltype" for T ∈ [ Int, Float64, Any ]
    x = convert(Vector{T}, collect(1:5))
    itp = linear_interpolation(x, x.^2)
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
        etp = linear_interpolation(1:10, rand(10),
            extrapolation_bc=((RevBC, FwdBC),)
        )
        kiter = knots(etp)
        @test typeof(kiter) <: Interpolations.KnotIterator
        ktype = KnotIterator{Int, Tuple{typeof(RevBC), typeof(FwdBC)}}
        @test typeof(kiter) == ktype

        @test_knots kiter KnotIterator knot_ref(1:10, FwdBC)

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

    @test_knots kiter KnotIterator 1:5

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

    @test_knots kiter KnotIterator 1:3 1:3

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
    @test_knots kiter KnotIterator 1:5

    # Non-repeating Knots -> Should iterate over once and be done
    @testset "extrapolation - $etp" for etp ∈ [ Throw(), Flat(), Line()]
        extrp = extrapolate(itp, etp)
        kiter = knots(extrp)
        @test_knots kiter KnotIterator knot_ref(1:5, etp)
    end

    # Repeating Knots -> Should Iterate indefinitely
    @testset "extrapolation - Periodic" begin
        extrp = extrapolate(itp, Periodic())
        k = knots(extrp)
        @test_knots k KnotIterator knot_ref(1:5, Periodic())

        # Check that knots maps to the correct inbound knot
        k10 = Iterators.take(k, 10) |> collect
        @test extrp.(k10) ≈ vcat(x[1:end-1].^2, x[1:end-1].^2, x[1:2].^2)
    end
    @testset "extrapolation - Reflect" begin
        extrp = extrapolate(itp, Reflect())
        k = knots(extrp)
        @test_knots k KnotIterator knot_ref(1:5, Reflect())


        # Check that knots maps to the correct inbound knot
        k10 = Iterators.take(k, 10) |> collect
        @test k10 ≈ collect(1:10)
        @test extrp.(k10) ≈ vcat(x.^2, reverse(x[1:end - 1].^2), x[2].^2)
    end
end

# Unit tests for iteration on an uneven grid - Periodic
# Knots should repeat indefinitely with the first and last knots being co-located
@testset "iterate - uneven - Periodic" begin
    x = [1.0, 1.3, 2.4, 3.2, 4.0]
    etp = linear_interpolation(x, x.^2, extrapolation_bc=Periodic())
    kiter = knots(etp)
    @test_knots kiter KnotIterator knot_ref(x, Periodic())

    # Check that knots maps to the correct inbound knot
    k = Iterators.take(kiter, 10) |> collect
    @test k == [1.0, 1.3, 2.4, 3.2, 4.0, 4.3, 5.4, 6.2, 7.0, 7.3]
    @test etp.(k) ≈ vcat(x[1:end-1], x[1:end-1], x[1:2]).^2
end

# Unit tests for iteration on an uneven grid - Reflect
@testset "iterate - uneven - Reflect" begin
    x = [1.0, 1.3, 2.4, 3.2, 4.0]
    etp = linear_interpolation(x, x.^2, extrapolation_bc=Reflect())
    kiter = knots(etp)
    @test_knots kiter KnotIterator knot_ref(x, Reflect())

    # Check that knots maps to the correct inbound knot
    k = Iterators.take(kiter, 10) |> collect
    @test k == [1.0, 1.3, 2.4, 3.2, 4.0, 4.8, 5.6, 6.7, 7.0, 7.3]
    @test etp.(k) ≈ vcat(x, reverse(x[1:end - 1]), x[2]).^2
end

# Unit tests for 2D iteration with directional boundary conditions that are
# bounded (ie. knots do not repeat indefinitely)
@testset "2D - iteration - bounded - $bc" for bc ∈ [Line(), (Throw(), Line())]
    etp = linear_interpolation(([1, 2, 3], [1, 2, 3]), rand(3, 3);
        extrapolation_bc=(Throw(), bc)
    )
    kiter = knots(etp)
    @test_knots kiter KnotIterator 1:3 1:3
end

# Unit tests for 2D iteration with directional boundary conditions that are
# unbounded (ie. knots do repeat indefinitely)
@testset "2D - iteration - Unbounded - $bc" for bc ∈ [Periodic(), (Throw(), Periodic())]
    etp = linear_interpolation(([1, 2, 3], [1, 2, 3]), rand(3, 3);
        extrapolation_bc=(Line(), bc)
    )
    kiter = knots(etp)
    @test_knots kiter KnotIterator 1:3 Iterators.countfrom(1)
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

@testset "knotsbetween - extrapolate - 1D - $etp" for etp ∈ [Throw(), Flat(), Line()]
    x = [1.0, 1.5, 1.75, 2.0]
    etp = linear_interpolation(x, x.^2, extrapolation_bc=etp)
    @testset "start and stop" begin
        krange = knotsbetween(etp; start=0.0, stop=4.2)
        @test_knots krange KnotRange [1.0, 1.5, 1.75, 2.0]
    end
    @testset "start" begin
        krange = knotsbetween(etp; start=1.2)
        @test_knots krange KnotRange [1.5, 1.75, 2.0]
    end
    @testset "stop" begin
        krange = knotsbetween(etp; start=1.2)
        @test_knots krange KnotRange [1.5, 1.75, 2.0]
    end
    @testset "empty iterator" begin
        krange = knotsbetween(etp; start=3)
        @test_knots krange KnotRange Float64[]

        krange = knotsbetween(etp; stop=1.0)
        @test_knots krange KnotRange Float64[]
    end
end

macro test_knot_idx(f, iter, start, state, knot)
    quote
        # Wrap in @testset to get debug info on failure
        local $iter = $(esc(iter))
        @testset "$($f)($($iter)), $($start)) -> $($knot), $($state)" begin
            # Check that state does correspond to knot
            knot = first(iterate($iter, $state))
            @test knot == $knot

            # Check that f(start) -> state
            state = $f($iter, $start)
            @test state == $state
        end
    end
end

@testset "_knot_start/stop - Periodic" begin
    x = [1.0, 1.5, 1.75, 2.0]
    etp = linear_interpolation(x, x.^2, extrapolation_bc=Periodic())
    krange = knotsbetween(etp; stop = 10)
    using Interpolations: nextknotidx, priorknotidx

    @test_knot_idx(nextknotidx, krange, -2.1, -8, -2.0)
    @test_knot_idx(nextknotidx, krange, 0.4, -1, 0.5)
    @test_knot_idx(nextknotidx, krange, 0.7, 0, 0.75)
    @test_knot_idx(nextknotidx, krange, 0.8, 1, 1.0)
    @test_knot_idx(nextknotidx, krange, 1.0, 2, 1.5)
    @test_knot_idx(nextknotidx, krange, 1.5, 3, 1.75)
    @test_knot_idx(nextknotidx, krange, 2.3, 5, 2.5)

    @test_knot_idx(priorknotidx, krange, 0.4, -2, 0.0)
    @test_knot_idx(priorknotidx, krange, 0.7, -1, 0.5)
    @test_knot_idx(priorknotidx, krange, 0.8, 0, 0.75)
    @test_knot_idx(priorknotidx, krange, 1.0, 0, 0.75)
    @test_knot_idx(priorknotidx, krange, 1.5, 1, 1.0)
    @test_knot_idx(priorknotidx, krange, 2.0, 3, 1.75)
    @test_knot_idx(priorknotidx, krange, 2.3, 4, 2.0)
end

@testset "_knot_start/stop - Reflect" begin
    x = [1.0, 1.5, 1.75, 2.0]
    etp = linear_interpolation(x, x.^2, extrapolation_bc=Reflect())
    krange = knotsbetween(etp; stop = 10)
    using Interpolations: nextknotidx, priorknotidx

    @test_knot_idx(nextknotidx, krange, -2.1, -8, -2.0)
    @test_knot_idx(nextknotidx, krange, 0.2, -1, 0.25)
    @test_knot_idx(nextknotidx, krange, 1.6, 3, 1.75)
    @test_knot_idx(nextknotidx, krange, 3.75, 10, 4.0)
    @test_knot_idx(nextknotidx, krange, 5.4, 14, 5.5)

    @test_knot_idx(priorknotidx, krange, 0.4, -1, 0.25)
    @test_knot_idx(priorknotidx, krange, 0.7, 0, 0.5)
    @test_knot_idx(priorknotidx, krange, 4.3, 11, 4.25)
    @test_knot_idx(priorknotidx, krange, 6.3, 17, 6.25)
    @test_knot_idx(priorknotidx, krange, 1.5, 1, 1.0)
    @test_knot_idx(priorknotidx, krange, 8.2, 22, 8.0)
    @test_knot_idx(priorknotidx, krange, 10.0, 27, 9.75)
end

@testset "knotsbetween - extrapolate - 1D - Periodic" begin
    x = [1.0, 1.5, 1.75, 2.0]
    etp = linear_interpolation(x, x.^2, extrapolation_bc=Periodic())
    @testset "start and stop" begin
        krange = knotsbetween(etp; start=0.0, stop=4.2)
        expect = knot_ref([0.5, 0.75, 1.0, 1.5], Periodic(); stop=4.2)
        @test_knots krange KnotRange expect
        @test all(0.0 .< collect(krange) .< 4.2)
    end
    @testset "start" begin
        krange = knotsbetween(etp; start=3.2)
        expect = knot_ref([3.5, 3.75, 4.0, 4.5], Periodic())
        @test_knots krange KnotRange expect
    end
    @testset "stop" begin
        krange = knotsbetween(etp; stop=4.8)
        expect = knot_ref([1.0, 1.5, 1.75, 2.0], Periodic(); stop=4.8)
        @test_knots krange KnotRange expect
        @test all(collect(krange) .< 4.8)
    end
    @testset "empty iterator" begin
        krange = knotsbetween(etp; start=1, stop=1)
        @test_knots krange KnotRange Float64[]
    end
end

@testset "knotsbetween - extrapolate - 1D - Reflect" begin
    x = [1.0, 1.5, 1.75, 2.0]
    etp = linear_interpolation(x, x.^2, extrapolation_bc=Reflect())
    @testset "start and stop" begin
        krange = knotsbetween(etp; start=0.0, stop=4.2)
        expect = knot_ref([-1.0, -0.5, -0.25, 0.0], Reflect(); start=0.0, stop=4.2)
        @test_knots krange KnotRange expect
        @test all(0.0 .< collect(krange) .< 4.2)
    end
    @testset "start" begin
        krange = knotsbetween(etp; start=3.2)
        expect = knot_ref([3.0, 3.5, 3.75, 4.0], Reflect(); start=3.2)
        @test_knots krange KnotRange expect
    end
    @testset "stop" begin
        krange = knotsbetween(etp; stop=4.8)
        expect = knot_ref([1.0, 1.5, 1.75, 2.0], Reflect(); stop=4.8)
        @test_knots krange KnotRange expect
        @test all(collect(krange) .< 4.8)
    end
    @testset "empty iterator" begin
        krange = knotsbetween(etp; start=1, stop=1)
        @test_knots krange KnotRange Float64[]
    end
end

@testset "knotsbetween - empty iterator" begin
    x = [1.0, 1.5, 1.75, 2.0]
    etp = linear_interpolation(x, x.^2)
    @testset "start is out of bounds" begin
        krange = knotsbetween(etp; start=2.0)
        @test_knots krange KnotRange Float64[]
    end
    @testset "start and stop and out of bounds" begin
        krange = knotsbetween(etp; start=-10, stop=0)
        @test_knots krange KnotRange Float64[]
    end
    @testset "stop is less than start" begin
        krange = knotsbetween(etp; start=1.8, stop=1.1)
        @test_knots krange KnotRange Float64[]
    end
end

@testset "knotsbetween - missing start and stop" begin
    @testset "1D Case" begin
        x = [1.0, 1.5, 1.75, 2.0]
        etp = linear_interpolation(x, x.^2)
        @test_throws ArgumentError knotsbetween(etp)
        @test_throws ArgumentError knotsbetween(etp; start=nothing, stop=nothing)
    end
    @testset "2D case" begin
        etp = linear_interpolation((1:3, 1:3), rand(3,3))
        @test_throws ArgumentError knotsbetween(etp)
        @test_throws ArgumentError knotsbetween(etp; start=(nothing, 1), stop=(nothing, 2))
    end
end

@testset "KnotIterator - Out of Bound Knots" begin
    x = [1.0, 1.5, 1.75, 2.0]
    @testset "NonRepeating Knots - $etp}" for etp ∈ [Flat(), Line(), Throw()]
        etp = linear_interpolation(x, x.^2, extrapolation_bc=etp)
        kiter = knots(etp)

        # Check that Bounds errors are thrown
        @test_throws BoundsError kiter[0]
        @test kiter[1] == x[1]
        @test kiter[kiter.nknots] == x[end]
        @test_throws BoundsError kiter[kiter.nknots+1]
    end
    @testset "Repeating Knots - $etp" for etp ∈ [Periodic(), Reflect()]
        etp = linear_interpolation(x, x.^2, extrapolation_bc=etp)
        kiter = knots(etp)

        # Check that Bounds errors are not thrown
        @test_nowarn kiter[0]
        @test kiter[1] == x[1]
        @test kiter[kiter.nknots] == x[end]
        @test_nowarn kiter[kiter.nknots+1]
    end
end
