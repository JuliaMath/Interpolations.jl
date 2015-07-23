module Issue34

using Interpolations, Base.Test

A = rand(1:20, 100, 100)

# In #34, this incantation throws
itp = interpolate(A, BSpline(Quadratic(Flat)), OnCell)
# Sanity check that not only don't throw, but actually interpolate
for i in 1:size(A,1), j in 1:size(A,2)
    @test_approx_eq itp[i,j] A[i,j]
end

end
