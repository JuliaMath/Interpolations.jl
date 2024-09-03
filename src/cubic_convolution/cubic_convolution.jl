export CubicConvolutionalInterpolation, ConvolutionMethod

# for type stability of specialized coefficient generation
const Knots1D = Tuple{AbstractVector{T}} where T
const Knots2D = Tuple{AbstractVector{T}, AbstractVector{T}} where T
const Knots3D = Tuple{AbstractVector{T}, AbstractVector{T}, AbstractVector{T}} where T

struct ConvolutionMethod <: InterpolationType end

# for type stability of specialized interpolation functions
struct HigherDimension{N} end
HigherDimension(::Val{N}) where N = HigherDimension{N}()
struct OrderOfAccuracy{O} end
OrderOfAccuracy(::Val{O}) where O = OrderOfAccuracy{O}()
struct CubicConvolutionalKernel{O} end
CubicConvolutionalKernel(::Val{O}) where O = CubicConvolutionalKernel{O}()

struct CubicConvolutionalInterpolation{T,N,TCoefs<:AbstractArray,IT<:NTuple{N,ConvolutionMethod},Axs<:Tuple,KA,DT,O} <: AbstractInterpolation{T,N,IT}
    coefs::TCoefs
    knots::Axs
    it::IT
    h::NTuple{N,Float64}
    kernel::KA
    dimension::DT
    order::O
end

include("cubic_convolution_coefs.jl")
include("cubic_convolution_kernels.jl")
include("cubic_convolution_extrapolation.jl")
include("cubic_convolution_interpolation.jl")

function CubicConvolutionalInterpolation(knots::NTuple{N,AbstractVector}, vs::AbstractArray{T,N}; order::Int=4) where {T,N}
    
    if order == 3
        coefs = create_cubic_convolutional_coefs(knots, vs) # expand boundaries once
    elseif order == 4
        coefs_init = create_cubic_convolutional_coefs(knots, vs)
        coefs = create_cubic_convolutional_coefs(knots, coefs_init) # expand boundaries twice
    else
        throw(ArgumentError("order must be 3 or 4"))
    end

    h = map(k -> k[2] - k[1], knots)
    it = ntuple(_ -> ConvolutionMethod(), N)

    if order == 3
        knots_new = expand_knots(knots, h) # expand boundaries once
    elseif order == 4
        knots_init = expand_knots(knots, h)
        knots_new = expand_knots(knots_init, h) # expand boundaries twice
    end

    dimension = N <= 3 ? Val(N) : HigherDimension(Val(N))
    cubickernel = CubicConvolutionalKernel(Val(order))
    bigO = OrderOfAccuracy(Val(order))

    CubicConvolutionalInterpolation{T,N,typeof(coefs),typeof(it),typeof(knots_new),typeof(cubickernel),typeof(dimension),typeof(bigO)}(
        coefs, knots_new, it, h, cubickernel, dimension, bigO
    )
end

function expand_knots(knots::NTuple{N,AbstractVector}, h::NTuple{N,Real}) where N
    knots_new = ntuple(i -> knots[i][1]-h[i]:h[i]:knots[i][end]+h[i], N)
    return knots_new
end

Interpolations.getknots(itp::CubicConvolutionalInterpolation) = itp.knots
Base.axes(itp::CubicConvolutionalInterpolation) = axes(itp.coefs)
Base.size(itp::CubicConvolutionalInterpolation) = size(itp.coefs)
Interpolations.lbounds(itp::CubicConvolutionalInterpolation) = first.(itp.knots)
Interpolations.ubounds(itp::CubicConvolutionalInterpolation) = last.(itp.knots)
Interpolations.itpflag(::Type{<:CubicConvolutionalInterpolation{T,N,TCoefs,IT}}) where {T,N,TCoefs,IT} = IT()
Interpolations.coefficients(itp::CubicConvolutionalInterpolation) = itp.coefs
Base.length(itp::CubicConvolutionalInterpolation) = length(itp.coefs)
Base.iterate(itp::CubicConvolutionalInterpolation, state=1) = state > length(itp) ? nothing : (itp[state], state+1)
itpflag(itp::CubicConvolutionalInterpolation) = itp.it
lbound(ax::AbstractRange, itp::ConvolutionMethod) = lbound(ax, itp.order)

# accessing the coefficients
function Base.getindex(itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,CubicConvolutionalKernel{O},DT,OrderOfAccuracy{O}}, I::Vararg{Integer,N}) where {T,N,TCoefs,IT,Axs,O,DT}
    return itp.coefs[I...]
end
function (itp::CubicConvolutionalInterpolation{T,1,TCoefs,IT,Axs,CubicConvolutionalKernel{O},Val{1},OrderOfAccuracy{O}})(i::Integer) where {T,TCoefs,IT,Axs,O}
    return itp.coefs[i]
end
function (itp::CubicConvolutionalInterpolation{T,N,TCoefs,IT,Axs,CubicConvolutionalKernel{O},DT,OrderOfAccuracy{O}})(I::Vararg{Integer,N}) where {T,N,TCoefs,IT,Axs,O,DT}
    return itp.coefs[I...]
end
