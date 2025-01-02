export ConvolutionInterpolation, FastConvolutionInterpolation, ConvolutionMethod, convolution_interpolation

# for type stability of specialized interpolation functions
struct HigherDimension{N} end
HigherDimension(::Val{N}) where N = HigherDimension{N}()
struct ConvolutionKernel{DG} end
ConvolutionKernel(::Val{DG}) where {DG} = ConvolutionKernel{DG}()
struct GaussianConvolutionKernel{B} end
GaussianConvolutionKernel(::Val{B}) where B = GaussianConvolutionKernel{B}()

struct ConvolutionInterpolation{T,N,TCoefs<:AbstractArray,IT<:NTuple{N,ConvolutionMethod},
                                Axs<:Tuple,KA,DT,DG,EQ} <: AbstractConvolutionInterpolation{T,N,TCoefs,IT,Axs,KA,DT,DG,EQ}
    coefs::TCoefs
    knots::Axs
    it::IT
    h::NTuple{N,Float64}
    kernel::KA
    dimension::DT
    deg::DG
    eqs::EQ
end

struct FastConvolutionInterpolation{T,N,TCoefs<:AbstractArray,IT<:NTuple{N,ConvolutionMethod},
                                Axs<:Tuple,KA,DT,DG,EQ,PR,KP} <: AbstractConvolutionInterpolation{T,N,TCoefs,IT,Axs,KA,DT,DG,EQ}
    coefs::TCoefs
    knots::Axs
    it::IT
    h::NTuple{N,Float64}
    kernel::KA
    dimension::DT
    deg::DG
    eqs::EQ
    pre_range::PR
    kernel_pre::KP
end

include("convolution_extrapolation.jl")
include("convolution_coefs.jl")
include("convolution_fast_interpolation.jl")
include("convolution_kernel_interpolation.jl")
include("convolution_kernels.jl")

function ConvolutionInterpolation(knots::NTuple{N,AbstractVector}, vs::AbstractArray{T,N};
                                degree::Int=3, B=nothing) where {T,N}

    eqs = B === nothing ? get_equations_for_degree(degree) : 50
    h = map(k -> k[2] - k[1], knots)
    it = ntuple(_ -> ConvolutionMethod(), N)

    knots_new = expand_knots(knots, eqs-1) # expand boundaries
    coefs = create_convolutional_coefs(vs, h, eqs) # create boundaries
    kernel = B === nothing ? ConvolutionKernel(Val(degree)) : GaussianConvolutionKernel(Val(B))
    dimension = N <= 3 ? Val(N) : HigherDimension(Val(N))
    degree = Val(degree)

    ConvolutionInterpolation{T,N,typeof(coefs),typeof(it),typeof(knots_new),typeof(kernel),typeof(dimension),typeof(degree),typeof(eqs)}(
        coefs, knots_new, it, h, kernel, dimension, degree, eqs
    )
end

function FastConvolutionInterpolation(knots::NTuple{N,AbstractVector}, vs::AbstractArray{T,N};
                                degree::Int=3, precompute::Int=1000, B=nothing) where {T,N}
    
    eqs = B === nothing ? get_equations_for_degree(degree) : 50
    h = map(k -> k[2] - k[1], knots)
    it = ntuple(_ -> ConvolutionMethod(), N)

    knots_new = expand_knots(knots, eqs-1) # expand boundaries
    coefs = create_convolutional_coefs(vs, h, eqs) # create boundaries
    kernel = B === nothing ? ConvolutionKernel(Val(degree)) : GaussianConvolutionKernel(Val(B))
    dimension = N <= 3 ? Val(N) : HigherDimension(Val(N))
    pre_range = range(0.0, 1.0, length=precompute)
    kernel_pre = zeros(T, precompute, 2*eqs)
    for i = 1:2*eqs
        kernel_pre[:,i] .= kernel.(pre_range .- eqs .+ i .- 1)
    end
    degree = Val(degree)

    FastConvolutionInterpolation{T,N,typeof(coefs),typeof(it),typeof(knots_new),typeof(kernel),typeof(dimension),typeof(degree),typeof(eqs),typeof(pre_range),typeof(kernel_pre)}(
        coefs, knots_new, it, h, kernel, dimension, degree, eqs, pre_range, kernel_pre
    )
end

function extend_vector(x::AbstractVector, n_extra::Integer)
    step_start = x[2] - x[1]
    step_end = x[end] - x[end-1]
    
    start_extension = range(x[1] - n_extra * step_start, step=step_start, length=n_extra)
    end_extension = range(x[end] + step_end, step=step_end, length=n_extra)
    
    return vcat(start_extension, x, end_extension)
end

function expand_knots(knots::NTuple{N,AbstractVector}, n_extra::Integer) where N
    knots_new = ntuple(i -> extend_vector(knots[i], n_extra), N)
    return knots_new
end

getknots(itp::ConvolutionInterpolation) = itp.knots
Base.axes(itp::ConvolutionInterpolation) = axes(itp.coefs)
Base.size(itp::ConvolutionInterpolation) = size(itp.coefs)
lbounds(itp::ConvolutionInterpolation) = first.(itp.knots)
ubounds(itp::ConvolutionInterpolation) = last.(itp.knots)
itpflag(::Type{<:ConvolutionInterpolation{T,N,TCoefs,IT}}) where {T,N,TCoefs,IT} = IT()
coefficients(itp::ConvolutionInterpolation) = itp.coefs
Base.length(itp::ConvolutionInterpolation) = length(itp.coefs)
Base.iterate(itp::ConvolutionInterpolation, state=1) = state > length(itp) ? nothing : (itp[state], state+1)
itpflag(itp::ConvolutionInterpolation) = itp.it
getknots(itp::FastConvolutionInterpolation) = itp.knots
Base.axes(itp::FastConvolutionInterpolation) = axes(itp.coefs)
Base.size(itp::FastConvolutionInterpolation) = size(itp.coefs)
lbounds(itp::FastConvolutionInterpolation) = first.(itp.knots)
ubounds(itp::FastConvolutionInterpolation) = last.(itp.knots)
itpflag(::Type{<:FastConvolutionInterpolation{T,N,TCoefs,IT}}) where {T,N,TCoefs,IT} = IT()
coefficients(itp::FastConvolutionInterpolation) = itp.coefs
Base.length(itp::FastConvolutionInterpolation) = length(itp.coefs)
Base.iterate(itp::FastConvolutionInterpolation, state=1) = state > length(itp) ? nothing : (itp[state], state+1)
itpflag(itp::FastConvolutionInterpolation) = itp.it
lbound(ax::AbstractRange, itp::ConvolutionMethod) = lbound(ax, itpflag(itp))

# accessing the coefficients
function Base.getindex(itp::ConvolutionInterpolation{T,N,TCoefs,IT,Axs,ConvolutionKernel{DG},DT}, I::Vararg{Integer,N}) where {T,N,TCoefs,IT,Axs,DT,DG}
    return itp.coefs[I...]
end
function (itp::ConvolutionInterpolation{T,1,TCoefs,IT,Axs,ConvolutionKernel{DG},Val{1}})(i::Integer) where {T,TCoefs,IT,Axs,DG}
    return itp.coefs[i]
end
function (itp::ConvolutionInterpolation{T,N,TCoefs,IT,Axs,ConvolutionKernel{DG},DT})(I::Vararg{Integer,N}) where {T,N,TCoefs,IT,Axs,DT,DG}
    return itp.coefs[I...]
end
function Base.getindex(itp::FastConvolutionInterpolation{T,N,TCoefs,IT,Axs,ConvolutionKernel{DG},DT}, I::Vararg{Integer,N}) where {T,N,TCoefs,IT,Axs,DT,DG}
    return itp.coefs[I...]
end
function (itp::FastConvolutionInterpolation{T,1,TCoefs,IT,Axs,ConvolutionKernel{DG},Val{1}})(i::Integer) where {T,TCoefs,IT,Axs,DG}
    return itp.coefs[i]
end
function (itp::FastConvolutionInterpolation{T,N,TCoefs,IT,Axs,ConvolutionKernel{DG},DT})(I::Vararg{Integer,N}) where {T,N,TCoefs,IT,Axs,DT,DG}
    return itp.coefs[I...]
end