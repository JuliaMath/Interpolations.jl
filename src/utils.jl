@inline sqr(x) = x*x
@inline cub(x) = x*x*x

modrange(x, r::AbstractUnitRange) = mod(x-first(r), length(r)) + first(r)
