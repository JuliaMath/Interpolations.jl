export ScaledInterpolation, eachvalue

import Base: iterate

struct ScaledInterpolation{T,N,ITPT,IT,RT} <: AbstractInterpolationWrapper{T,N,ITPT,IT}
    itp::ITPT
    ranges::RT
end

Base.parent(A::ScaledInterpolation) = A.itp
count_interp_dims(::Type{<:ScaledInterpolation{T,N,ITPT}}, n) where {T,N,ITPT} = count_interp_dims(ITPT, n)

"""
`scale(itp, xs, ys, ...)` scales an existing interpolation object to allow for indexing using other coordinate axes than unit ranges, by wrapping the interpolation object and transforming the indices from the provided axes onto unit ranges upon indexing.

The parameters `xs` etc must be either ranges or linspaces, and there must be one coordinate range/linspace for each dimension of the interpolation object.

For every `NoInterp` dimension of the interpolation object, the range must be exactly `1:size(itp, d)`.
"""
function scale(itp::AbstractInterpolation{T,N,IT}, ranges::Vararg{AbstractRange,N}) where {T,N,IT}
    check_ranges(itpflag(itp), axes(itp), ranges)
    ScaledInterpolation{T,N,typeof(itp),IT,typeof(ranges)}(itp, ranges)
end

function check_ranges(flags, axs, ranges)
    check_range(getfirst(flags), axs[1], ranges[1])
    check_ranges(getrest(flags), Base.tail(axs), Base.tail(ranges))
end
check_ranges(::Any, ::Tuple{}, ::Tuple{}) = nothing

check_range(::NoInterp, ax, r) = ax == r || throw(ArgumentError("The range $r did not equal the corresponding axis of the interpolation object $ax"))
check_range(::Any, ax, r) = length(ax) == length(r) || throw(ArgumentError("The range $r is incommensurate with the corresponding axis $ax"))

size(sitp::ScaledInterpolation) = size(sitp.itp)
axes(sitp::ScaledInterpolation) = axes(sitp.itp)

lbounds(sitp::ScaledInterpolation) = _lbounds(sitp.ranges, itpflag(sitp.itp), gridflag(sitp.itp))
ubounds(sitp::ScaledInterpolation) = _ubounds(sitp.ranges, itpflag(sitp.itp), gridflag(sitp.itp))

lbound(ax, ::BSpline, ::OnCell) = first(ax) - boundstep(ax)
ubound(ax, ::BSpline, ::OnCell) = last(ax) + boundstep(ax)
lbound(ax, ::BSpline, ::OnGrid) = first(ax)
ubound(ax, ::BSpline, ::OnGrid) = last(ax)

boundstep(r::StepRange) = r.step / 2
boundstep(r::UnitRange) = 1//2

"""
Returns *half* the width of one step of the range.

This function is used to calculate the upper and lower bounds of `OnCell` interpolation objects.
""" boundstep

function (sitp::ScaledInterpolation{T,N})(xs::Vararg{Number,N}) where {T,N}
    xl = coordslookup(itpflag(sitp.itp), sitp.ranges, xs)
    sitp.itp(xl...)
end

(sitp::ScaledInterpolation{T,1}, x::Number, y::Int) where {T} = y == 1 ? sitp(x) : Base.throw_boundserror(sitp, (x, y))

@inline function coordslookup(flags, ranges, xs)
    item = coordlookup(getfirst(flags), ranges[1], xs[1])
    (item, coordslookup(getrest(flags), Base.tail(ranges), Base.tail(xs))...)
end
coordslookup(::Any, ::Tuple{}, ::Tuple{}) = ()

coordlookup(::NoInterp, r, i) = i
coordlookup(::Flag, r, x) = coordlookup(r, x)

coordlookup(r::UnitRange, x) = x - r.start + oneunit(eltype(r))
# coordlookup(i::Bool, r::AbstractRange, x) = i ? coordlookup(r, x) : convert(typeof(coordlookup(r,x)), x)
coordlookup(r::StepRange, x) = (x - r.start) / r.step + oneunit(eltype(r))

coordlookup(r::StepRangeLen, x) = (x - first(r)) / step(r) + oneunit(eltype(r))
boundstep(r::StepRangeLen) = 0.5*step(r)
rescale_gradient(r::StepRangeLen, g) = g / step(r)

basetype(::Type{ScaledInterpolation{T,N,ITPT,IT,RT}}) where {T,N,ITPT,IT,RT} = ITPT
basetype(sitp::ScaledInterpolation) = basetype(typeof(sitp))


function gradient(sitp::ScaledInterpolation{T,N}, xs::Vararg{Number,N}) where {T,N}
    xl = coordslookup(itpflag(sitp.itp), sitp.ranges, xs)
    g = gradient(sitp.itp, xl...)
    SVector(rescale_gradient_components(itpflag(sitp.itp), sitp.ranges, Tuple(g)))
end

function rescale_gradient_components(flags, ranges, g)
    if getfirst(flags) isa NoInterp
        return rescale_gradient_components(getrest(flags), Base.tail(ranges), g)    # don't consume a coordinate of g
    else
        item = rescale_gradient(ranges[1], g[1])
        return (item, rescale_gradient_components(getrest(flags), Base.tail(ranges), Base.tail(g))...)
    end
end
rescale_gradient_components(flags, ::Tuple{}, ::Tuple{}) = ()


# # @eval uglyness required for disambiguation with method in b-splies/indexing.jl
# # also, GT is only specified to avoid disambiguation warnings on julia 0.4
# gradient(sitp::ScaledInterpolation{T,N,ITPT,IT}, xs::Real...) where {T,N,ITPT,IT<:DimSpec{InterpolationType}<:DimSpec{GridType}} =
#         gradient!(Array{T}(undef, count_interp_dims(IT,N)), sitp, xs...)
# gradient(sitp::ScaledInterpolation{T,N,ITPT,IT}, xs...) where {T,N,ITPT,IT<:DimSpec{InterpolationType}<:DimSpec{GridType}} =
#         gradient!(Array{T}(undef, count_interp_dims(IT,N)), sitp, xs...)
# @generated function gradient!(g, sitp::ScaledInterpolation{T,N,ITPT,IT}, xs::Number...) where {T,N,ITPT,IT}
#     ndims(g) == 1 || throw(DimensionMismatch("g must be a vector (but had $(ndims(g)) dimensions)"))
#     length(xs) == N || throw(DimensionMismatch("Must index into $N-dimensional scaled interpolation object with exactly $N indices (you used $(length(xs)))"))

#     interp_types = length(IT.parameters) == N ? IT.parameters : tuple([IT.parameters[1] for _ in 1:N]...)
#     interp_dimens = map(it -> interp_types[it] != NoInterp, 1:N)
#     interp_indices = map(i -> interp_dimens[i] ? :(coordlookup(sitp.ranges[$i], xs[$i])) : :(xs[$i]), 1:N)

#     quote
#         length(g) == $(count_interp_dims(IT, N)) || throw(ArgumentError(string("The length of the provided gradient vector (", length(g), ") did not match the number of interpolating dimensions (", $(count_interp_dims(IT, N)), ")")))
#         gradient!(g, sitp.itp, $(interp_indices...))
#         cntr = 0
#         for i = 1:N
#                 if $(interp_dimens)[i]
#                     cntr += 1
#                     g[cntr] = rescale_gradient(sitp.ranges[i], g[cntr])
#                 end
#         end
#         g
#     end
# end


rescale_gradient(r::StepRange, g) = g / r.step
rescale_gradient(r::UnitRange, g) = g

"""
`rescale_gradient(r::AbstractRange)`

Implements the chain rule dy/dx = dy/du * du/dx for use when calculating gradients with scaled interpolation objects.
""" rescale_gradient


# ### Iteration
# mutable struct ScaledIterator{CR<:CartesianIndices,SITPT,X1,Deg,T}
#     rng::CR
#     sitp::SITPT
#     dx_1::X1
#     nremaining::Int
#     fx_1::X1
#     itp_tail::NTuple{Deg,T}
# end

# nelements(::Union{Type{NoInterp},Type{Constant}}) = 1
# nelements(::Type{Linear}) = 2
# nelements(::Type{Q}) where {Q<:Quadratic} = 3

# eachvalue_zero(::Type{R}, ::Type{BT}) where {R,BT<:Union{Type{NoInterp},Type{Constant}}} =
#     (zero(R),)
# eachvalue_zero(::Type{R}, ::Type{Linear}) where {R} = (zero(R),zero(R))
# eachvalue_zero(::Type{R}, ::Type{Q}) where {R,Q<:Quadratic} = (zero(R),zero(R),zero(R))

# """
# `eachvalue(sitp)` constructs an iterator for efficiently visiting each
# grid point of a ScaledInterpolation object in which a small grid is
# being "scaled up" to a larger one.  For example, suppose you have a
# core `BSpline` object defined on a 5x7x4 grid, and you are scaling it
# to a 100x120x20 grid (via `linspace(1,5,100), linspace(1,7,120),
# linspace(1,4,20)`).  You can perform interpolation at each of these
# grid points via

# ```
#     function foo!(dest, sitp)
#         i = 0
#         for s in eachvalue(sitp)
#             dest[i+=1] = s
#         end
#         dest
#     end
# ```

# which should be more efficient than

# ```
#     function bar!(dest, sitp)
#         for I in CartesianIndices(size(dest))
#             dest[I] = sitp[I]
#         end
#         dest
#     end
# ```
# """
# function eachvalue(sitp::ScaledInterpolation{T,N}) where {T,N}
#     ITPT = basetype(sitp)
#     IT = itptype(ITPT)
#     R = getindex_return_type(ITPT, Int)
#     BT = bsplinetype(iextract(IT, 1))
#     itp_tail = eachvalue_zero(R, BT)
#     dx_1 = coordlookup(sitp.ranges[1], 2) - coordlookup(sitp.ranges[1], 1)
#     ScaledIterator(CartesianIndices(ssize(sitp)), sitp, dx_1, 0, zero(dx_1), itp_tail)
# end

# function index_gen1(::Union{Type{NoInterp}, Type{BSpline{Constant}}})
#     quote
#         value = iter.itp_tail[1]
#     end
# end

# function index_gen1(::Type{BSpline{Linear}})
#     quote
#         p = iter.itp_tail
#         value = c_1*p[1] + cp_1*p[2]
#     end
# end

# function index_gen1(::Type{BSpline{Q}}) where Q<:Quadratic
#     quote
#         p = iter.itp_tail
#         value = cm_1*p[1] + c_1*p[2] + cp_1*p[3]
#     end
# end
# function index_gen_tail(B::Union{Type{NoInterp}, Type{BSpline{Constant}}}, ::Type{IT}, N) where IT
#     [index_gen(B, IT, N, 0)]
# end

# function index_gen_tail(::Type{BSpline{Linear}}, ::Type{IT}, N) where IT
#     [index_gen(BS1, IT, N, i) for i = 0:1]
# end

# function index_gen_tail(::Type{BSpline{Q}}, ::Type{IT}, N) where {IT,Q<:Quadratic}
#     [index_gen(BSpline{Q}, IT, N, i) for i = -1:1]
# end
# function nremaining_gen(::Union{Type{BSpline{Constant}}, Type{BSpline{Q}}}) where Q<:Quadratic
#     quote
#         EPS = 0.001*iter.dx_1
#         floor(Int, iter.dx_1 >= 0 ?
#               (min(length(range1)+EPS, round(Int,x_1) + 0.5) - x_1)/iter.dx_1 :
#               (max(1-EPS, round(Int,x_1) - 0.5) - x_1)/iter.dx_1)
#     end
# end

# function nremaining_gen(::Type{BSpline{Linear}})
#     quote
#         EPS = 0.001*iter.dx_1
#         floor(Int, iter.dx_1 >= 0 ?
#               (min(length(range1)+EPS, floor(Int,x_1) + 1) - x_1)/iter.dx_1 :
#               (max(1-EPS, floor(Int,x_1)) - x_1)/iter.dx_1)
#     end
# end
# function next_gen(::Type{ScaledIterator{CR,SITPT,X1,Deg,T}}) where {CR,SITPT,X1,Deg,T}
#     N = ndims(CR)
#     ITPT = basetype(SITPT)
#     IT = itptype(ITPT)
#     BS1 = iextract(IT, 1)
#     BS1 == NoInterp && error("eachvalue is not implemented (and does not make sense) for NoInterp along the first dimension")
#     pad = padding(ITPT)
#     x_syms = [Symbol("x_", i) for i = 1:N]
#     interp_index(IT, i) = iextract(IT, i) != NoInterp ?
#         :($(x_syms[i]) = coordlookup(sitp.ranges[$i], state[$i])) :
#         :($(x_syms[i]) = state[$i])
#     # Calculations for the first dimension
#     interp_index1 = interp_index(IT, 1)
#     indices1 = define_indices_d(BS1, 1, padextract(pad, 1))
#     coefexprs1 = coefficients(BS1, N, 1)
#     nremaining_expr = nremaining_gen(BS1)
#     # Calculations for the rest of the dimensions
#     interp_indices_tail = map(i -> interp_index(IT, i), 2:N)
#     indices_tail = [define_indices_d(iextract(IT, i), i, padextract(pad, i)) for i = 2:N]
#     coefexprs_tail = [coefficients(iextract(IT, i), N, i) for i = 2:N]
#     value_exprs_tail = index_gen_tail(BS1, IT, N)
#     quote
#         sitp = iter.sitp
#         itp = sitp.itp
#         inds_itp = axes(itp)
#         if iter.nremaining > 0
#             iter.nremaining -= 1
#             iter.fx_1 += iter.dx_1
#         else
#             range1 = sitp.ranges[1]
#             $interp_index1
#             $indices1
#             iter.nremaining = $nremaining_expr
#             iter.fx_1 = fx_1
#             $(interp_indices_tail...)
#             $(indices_tail...)
#             $(coefexprs_tail...)
#             @inbounds iter.itp_tail = ($(value_exprs_tail...),)
#         end
#         fx_1 = iter.fx_1
#         $coefexprs1
#         $(index_gen1(BS1))
#     end
# end

# @generated function iterate(iter::ScaledIterator{CR,ITPT}, state::Union{Nothing,CartesianIndex{N}} = nothing) where {CR,ITPT,N}
#     value_expr = next_gen(iter)
#     quote
#         rng_next = state ≡ nothing ? iterate(iter.rng) : iterate(iter.rng, state)
#         rng_next ≡ nothing && return nothing
#         state = rng_next[2]
#         $value_expr
#         (value, state)
#     end
# end

# ssize(sitp::ScaledInterpolation{T,N}) where {T,N} = map(r->round(Int, last(r)-first(r)+1), sitp.ranges)::NTuple{N,Int}
