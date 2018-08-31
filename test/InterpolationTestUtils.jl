module InterpolationTestUtils

using Test, Interpolations, ForwardDiff, StaticArrays
using Interpolations: degree, itpflag, bounds, lbounds, ubounds

export check_axes, check_inbounds_values, check_oob, can_eval_near_boundaries,
       check_gradient, check_hessian
export MyPair

const failstore = Ref{Any}(nothing)   # stash the inputs to failing tests here

## Property accessors
coefs(itp) = itp.coefs

boundaryconditions(itp::AbstractInterpolation) = boundaryconditions(itpflag(itp))
boundaryconditions(bs::BSpline) = degree(bs)
boundaryconditions(ni::NoInterp) = ni

ndaccessor(x, d) = x
ndaccessor(x::Tuple, d) = x[d]

haspadding(itp) = haspadding(boundaryconditions(itp))
haspadding(::Union{Constant,Linear,NoInterp}) = false
haspadding(::Quadratic{BC}) where BC = haspadding(BC())
haspadding(::Cubic{BC})     where BC = haspadding(BC())
haspadding(::BC) where BC<:Interpolations.BoundaryCondition = !(BC <: Union{Periodic,InPlace,InPlaceQ})
haspadding(bcs::Tuple) = map(haspadding, bcs)
haspadding(::NoInterp) = false

getindexib(itp, i...) = @inbounds itp[i...]
callib(itp, i...)     = @inbounds itp(i...)

⊂(r1::AbstractRange, r2::AbstractRange) = first(r2) < first(r1) < last(r2) && first(r2) < last(r1) < last(r2)

function check_axes(itp, A, isinplace=false)
    @test ndims(itp) == ndims(A)
    axsi, axsA = @inferred(axes(itp)), axes(A)
    szi, szA = @inferred(size(itp)), size(A)
    haspad = haspadding(itp)
    for d = 1:ndims(A)
        if isinplace && ndaccessor(haspad, d)
            @test axsi[d] != axsA[d] && axsi[d] ⊂ axsA[d]
            @test szi[d] < szA[d]
        else
            @test axsi[d] == axsA[d]
            @test szi[d] == szA[d]
        end
    end
    nothing
end

function check_inbounds_values(itp, A)
    i = first(eachindex(itp))
    @test A[i] ≈ @inferred(itp[i]) == @inferred(itp[Tuple(i)...]) ≈ @inferred(itp(i)) ≈ @inferred(itp(float.(Tuple(i))...))
    for i in eachindex(itp)
        @test A[i] ≈ itp[i] == itp[Tuple(i)...] ≈ itp(i) ≈ itp(float.(Tuple(i))...)
    end
    cb = Base.JLOptions().check_bounds == 1
    if ndims(itp) == 1
        for i in eachindex(itp)
            @test itp[i,1] ≈ A[i]   # used in the AbstractArray display infrastructure
            @test_throws BoundsError itp[i,2]
            # if cb
            #     @test_throws BoundsError getindexib(itp, i, 2)
            #     @test_throws BoundsError callib(itp, i, 2)
            # else
            #     @test getindexib(itp, i, 2) ≈ A[i]
            #     @test callib(itp, i, 2) ≈ A[i]
            # end
        end
    end
    nothing
end

function check_oob(itp)
    widen(r) = first(r)-1:last(r)+1
    indsi = axes(itp)
    indsci = CartesianIndices(indsi)
    for i in CartesianIndices(widen.(indsi))
        i ∈ indsci && continue
        @test_throws BoundsError itp[i]
    end
    nothing
end

function can_eval_near_boundaries(itp::AbstractInterpolation{T,1}) where T
    l, u = bounds(itp, 1)
    @test isfinite(itp(l+0.1))
    @test isfinite(itp(u-0.1))
    @test_throws BoundsError itp(l-0.1)
    @test_throws BoundsError itp(u+0.1)
end

function can_eval_near_boundaries(itp::AbstractInterpolation)
    l, u = Float64.(lbounds(itp)), Float64.(ubounds(itp))
    for d = 1:ndims(itp)
        nearl = substitute(l, d, l[d]+0.1)
        # @show summary(itp) nearl
        @test isfinite(itp(nearl...))
        outl = substitute(l, d, l[d]-0.1)
        @test_throws BoundsError itp(outl...)
        nearu = substitute(u, d, u[d]-0.1)
        # @show nearu
        @test isfinite(itp(nearu...))
        outu = substitute(u, d, u[d]+0.1)
        @test_throws BoundsError itp(outu...)
    end
end

function substitute(default::NTuple{N,T}, d::Integer, val::T) where {T,N}
    ntuple(i->ifelse(i==d, val, default[i]), Val(N))
end


# Generate a grid of points [1.0, 1.3333, 1.6667, 2.0, 2.3333, ...] along each coordinate
thirds(axs) = Iterators.product(_thirds(axs...)...)

_thirds(a, axs...) =
    (sort(Float64[a; (first(a):last(a)-1) .+ 1/3; (first(a)+1:last(a)) .- 1/3]), _thirds(axs...)...)
_thirds() = ()

function check_gradient(itp::AbstractInterpolation, gtmp)
    val(x) = itp(Tuple(x)...)
    g!(gstore, x) = ForwardDiff.gradient!(gstore, val, x)
    gtmp2 = similar(gtmp)
    i = first(thirds(axes(itp)))
    @inferred(Interpolations.gradient(itp, i...))
    for i in thirds(axes(itp))
        @test Interpolations.gradient(itp, i...) ≈ g!(gtmp, SVector(i))
        @test Interpolations.gradient!(gtmp2, itp, i...) ≈ gtmp
    end
end

function check_hessian(itp::AbstractInterpolation, htmp)
    val(x) = itp(Tuple(x)...)
    h!(hstore, x) = ForwardDiff.hessian!(hstore, val, x)
    htmp2 = similar(htmp)
    i = first(thirds(axes(itp)))
    @inferred(Interpolations.hessian(itp, i...))
    for i in thirds(axes(itp))
        @test Interpolations.hessian(itp, i...) ≈ h!(htmp, SVector(i))
        @test Interpolations.hessian!(htmp2, itp, i...) ≈ htmp
    end
end

## A type used for multi-valued tests
import Base: +, -, *, /, ≈

struct MyPair{T}
    first::T
    second::T
end

# Here's the interface your type must define
(+)(p1::MyPair, p2::MyPair) = MyPair(p1.first+p2.first, p1.second+p2.second)
(-)(p1::MyPair, p2::MyPair) = MyPair(p1.first-p2.first, p1.second-p2.second)
(*)(n::Number, p::MyPair) = MyPair(n*p.first, n*p.second)
(*)(p::MyPair, n::Number) = n*p
(/)(p::MyPair, n::Number) = MyPair(p.first/n, p.second/n)
Base.zero(::Type{MyPair{T}}) where {T} = MyPair(zero(T),zero(T))
Base.promote_rule(::Type{MyPair{T1}}, ::Type{T2}) where {T1,T2<:Number} = MyPair{promote_type(T1,T2)}
≈(p1::MyPair, p2::MyPair) = (p1.first ≈ p2.first) & (p1.second ≈ p2.second)

end
