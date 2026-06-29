module InterpolationsSciMLBaseExt

using Interpolations: AbstractInterpolation
using SciMLBase: SciMLBase

# desired behavior is to warn if the contained abstractarray would warn
# assumes that the interpolation is not empty.
function SciMLBase.should_warn_paramtype(i::AbstractInterpolation{<:AbstractArray})
    return SciMLBase.should_warn_paramtype(first(i))
end

end
