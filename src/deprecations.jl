# deprecate getindex for non-integer numeric indices
@deprecate getindex(itp::AbstractInterpolation{T,N}, i::Vararg{Number,N}) where {T,N} itp(i...)
@deprecate getindex(itp::AbstractInterpolation{T,N}, i::Vararg{ExpandedIndexTypes,N}) where {T,N} itp(i...)

for T in (:Throw, :Flat, :Line, :Free, :Periodic, :Reflect, :InPlace, :InPlaceQ)
    @eval begin
        # Support changing the gridtype for an instance
        (::$T)(gt::GridType) = $T(gt)
        $T{GT}() where GT<:GridType = $T(GT())
    end
end

@deprecate interpolate(A::AbstractArray, ::NoInterp, ::GridType) interpolate(A, NoInterp())
function interpolate(::Type{TWeights}, ::Type{TC}, A, it::IT, gt::GT) where {TWeights,TC,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    bcs = create_bcs(it, gt)
    Base.depwarn("interpolate($TWeights, $TC, A, $it, $gt) is deprecated, use interpolate($TWeights, $TC, A, $bcs)", :interpolate)
    interpolate(TWeights, TC, A, bcs)
end
function interpolate(A::AbstractArray, it::IT, gt::GT) where {IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    bcs = create_bcs(it, gt)
    Base.depwarn("interpolate(A, $it, $gt) is deprecated, use interpolate(A, $bcs)", :interpolate)
    interpolate(A, bcs)
end

function interpolate!(::Type{TWeights}, A, it::IT, gt::GT) where {TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    bcs = create_bcs(it, gt)
    Base.depwarn("interpolate!($TWeights, A, $it, $gt) is deprecated, use interpolate!($TWeights, A, $bcs)", :interpolate)
    interpolate!(TWeights, A, bcs)
end
function interpolate!(A, it::IT, gt::GT) where {TWeights,IT<:DimSpec{BSpline},GT<:DimSpec{GridType}}
    bcs = create_bcs(it, gt)
    Base.depwarn("interpolate!(A, $it, $gt) is deprecated, use interpolate!(A, $bcs)", :interpolate)
    interpolate!(A, bcs)
end

# extrapolate(A, Linear()) should probably have been extrapolate(A, Line())
# (Line<:BoundaryCondition but Linear<:Degree)
# @deprecate extrapolate(itp::AbstractInterpolation{T,N,IT}, ::Linear) where {T,N,IT} extrapolate(itp, Line())

const OldExtrapDimSpec = Union{Flag,Tuple{Vararg{Union{Flag,NTuple{2,Flag}}}}}
function extrapolate(itp::AbstractInterpolation{T,N,IT}, etpflag::OldExtrapDimSpec) where {T,N,IT}
    replacement = replace_linear_line(etpflag)
    io = IOBuffer()
    show(io, etpflag)
    etpstring = String(take!(io))
    show(io, replacement)
    repstring = String(take!(io))
    Base.depwarn("extrapolate(itp, $etpstring) is deprecated, use extrapolate(itp, $repstring) instead", :extrapolate)
    extrapolate(itp, replacement)
end

create_bcs(it::BSpline, gt::GridType) = BSpline(create_bcs(degree(it), gt))
create_bcs(it::NoInterp, gt::GridType) = it
create_bcs(it::Constant, gt::GridType) = it
create_bcs(it::Linear, gt::GridType) = it
create_bcs(it::Quadratic, gt::GridType) = it(gt)
create_bcs(it::Cubic, gt::GridType) = it(gt)
create_bcs(it::Tuple, gt::Tuple) = map((i,g)->i(g), it, gt)
create_bcs(it::Tuple, gt::GridType) = map(t->create_bcs(t, gt), it)
create_bcs(it::Flag, gt::Tuple) = map(t->it(t), gt)

replace_linear_line(::Linear) = Line()
replace_linear_line(bc::BoundaryCondition) = bc
replace_linear_line(etpflag::Tuple) = replace_linear_line.(etpflag)
