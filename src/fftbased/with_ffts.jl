import FFTW: irfft,rfft

struct FFT <: InterpolationType
    new_size::Number
end

function interpolate(sig::AbstractArray{<:Real,1},it::FFT)
    ## For real signal, so we can used the real fast fourier transform (rfft) to interpolate the signal
    length(sig) >= it.new_size && error("New size of signal must be superior or equal")
    n_orig = div(length(sig),2)
    return irfft( vcat(rfft(sig),zeros( div(it.new_size,2 ) - n_orig )  ), it.new_size)
end
