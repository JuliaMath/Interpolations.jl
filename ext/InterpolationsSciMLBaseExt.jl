module InterpolationsSciMLBaseExt

using Interpolations: AbstractInterpolation, coefficients
using SciMLBase: SciMLBase

# desired behavior is to warn if the contained array would warn
# so we just forward the coefficients of the interpolation to the actual warning check
function SciMLBase.should_warn_paramtype(i::AbstractInterpolation)
    return SciMLBase.should_warn_paramtype(coefficients(i))
end

end
