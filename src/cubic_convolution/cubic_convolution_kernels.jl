function (::CubicConvolutionalKernel{3})(s)
    s_abs = abs(s)
    if s_abs < 1.0
        return (3/2)*s_abs^3 - (5/2)*s_abs^2 + 1
    elseif s_abs < 2.0
        return -(1/2)*s_abs^3 + (5/2)*s_abs^2 - 4*s_abs + 2
    else
        return 0.0
    end
end

function (::CubicConvolutionalKernel{4})(s)
    s_abs = abs(s)
    if s_abs < 1.0
        return (4/3)*s_abs^3 - (7/3)*s_abs^2 + 1
    elseif s_abs < 2.0
        return -(7/12)*s_abs^3 + 3*s_abs^2 - (59/12)*s_abs + 15/6
    elseif s_abs < 3.0
        return (1/12)*s_abs^3 - (2/3)*s_abs^2 + (21/12)*s_abs - 3/2
    else
        return 0.0
    end
    
end