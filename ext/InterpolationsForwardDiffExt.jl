module InterpolationsForwardDiffExt

import Interpolations
using ForwardDiff

# this strips arbitrary layers of ForwardDiff.Dual, returning the innermost value
Interpolations.just_dual_value(x::ForwardDiff.Dual) = Interpolations.just_dual_value(ForwardDiff.value(x))

function Interpolations.maybe_clamp(::Interpolations.NeedsCheck, itp, xs::Tuple{Vararg{ForwardDiff.Dual}})
    xs_values = Interpolations.just_dual_value.(xs)
    clamped_vals = Interpolations.maybe_clamp(Interpolations.NeedsCheck(), itp, xs_values)
    apply_partials.(xs, clamped_vals)
end

# apply partials from arbitrarily nested ForwardDiff.Dual to a value
# used in maybe_clamp, above
function apply_partials(dual_number::D, x::Number) where D <: ForwardDiff.Dual
    ∂s = ForwardDiff.partials(dual_number)
    val = ForwardDiff.value(dual_number)
    D(apply_partials(val, x) , ∂s)
end
apply_partials(_::Number, x::Number) = x

end
