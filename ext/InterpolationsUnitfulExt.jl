module InterpolationsUnitfulExt

import Unitful, Interpolations

Interpolations.tweight(A::AbstractArray{T}) where T <: Unitful.Quantity = Interpolations.tweight(Unitful.ustrip(A))

end
