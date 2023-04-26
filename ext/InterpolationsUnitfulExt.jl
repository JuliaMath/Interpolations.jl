module InterpolationsUnitfulExt

if !isdefined(Base, :get_extension)
    import ..Unitful, ..Interpolations
else
    import Unitful, Interpolations
end

Interpolations.tweight(A::AbstractArray{T}) where T <: Unitful.Quantity = Interpolations.tweight(Unitful.ustrip(A))

end
