import .Unitful

tweight(A::AbstractArray{T}) where T <: Unitful.Quantity = tweight(Unitful.ustrip(A))