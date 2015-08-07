using Base.Cartesian

import Base.getindex

# Indexing at a point
@generated function getindex{T,N,TCoefs,IT<:DimSpec{Gridded},K,P}(itp::GriddedInterpolation{T,N,TCoefs,IT,K,P}, x::Number...)
    length(x) == N || error("Can only be called with $N indexes")
    meta = Expr(:meta, :inline)
    quote
        $meta
        @nexprs $N d->begin
            x_d = x[d]
            k_d = itp.knots[d]
            ix_d = searchsortedfirst(k_d, x_d, 1, length(k_d), Base.Order.ForwardOrdering()) - 1
        end
        $(define_indices(IT, N, P))
        $(coefficients(IT, N))
        @inbounds ret = $(index_gen(IT, N))
        ret
    end
end

# Because of the "vectorized" definition below, we need a definition for CartesianIndex
@generated function getindex{T,N}(itp::GriddedInterpolation{T,N}, index::CartesianIndex{N})
    args = [:(index[$d]) for d = 1:N]
    :(getindex(itp, $(args...)))
end

# Indexing with vector inputs. Here, it pays to pre-process the input indexes,
# because N*n is much smaller than n^N.
# TODO: special-case N=1, because there is no reason to separately cache the indexes.
@generated function getindex!{T,N,TCoefs,IT<:DimSpec{Gridded},K,P}(dest, itp::GriddedInterpolation{T,N,TCoefs,IT,K,P}, xv...)
    length(xv) == N || error("Can only be called with $N indexes")
    indexes_exprs = Expr[define_indices_d(iextract(IT, d), d, P) for d = 1:N]
    coefficient_exprs = Expr[coefficients(iextract(IT, d), N, d) for d = 1:N]
    # A manual @nloops (the interaction of d with the two exprs above is tricky...)
    ex = :(@inbounds @nref($N,dest,i) = $(index_gen(IT, N)))
    for d = 1:N
        isym, xsym, xvsym, ixsym, ixvsym = symbol("i_",d), symbol("x_",d), symbol("xv_",d), symbol("ix_",d), symbol("ixv_",d)
        ex = quote
            for $isym = 1:length($xvsym)
                $xsym  = $xvsym[$isym]
                $ixsym = $ixvsym[$isym]
                $(indexes_exprs[d])
                $(coefficient_exprs[d])
                $ex
            end
        end
    end
    quote
        @nexprs $N d->begin
            xv_d = xv[d]
            k_d = itp.knots[d]
            ixv_d = Array(Int, length(xv_d))  # ixv_d[i] is the smallest value such that k_d[ixv_d[i]] <= x_d[i]
            # If x_d is sorted and has quite a few entries, it's better to match
            # entries of x_d and k_d by iterating through them both in unison.
            l_d = length(k_d)   # FIXME: check l_d == 1 someday, see FIXME above
            # estimate the time required for searchsortedfirst vs. linear traversal
            den = 5*log(l_d) - 1   # 5 is arbitrary, for now (it's the coefficient of ssf compared to the while loop below)
            ascending = den*length(xv_d) > l_d  # if this is (or becomes) false, use searchsortedfirst
            i = 2  # this clamps ixv_d .>= 1
            knext = k_d[i]
            xjold = xv_d[1]
            for j = 1:length(xv_d)
                xj = xv_d[j]
                ascending = ascending & (xj >= xjold)
                if ascending
                    while i < length(k_d) && knext < xj
                        knext = k_d[i+=1]
                    end
                    ixv_d[j] = i-1
                    xjold = xj
                else
                    ixv_d[j] = searchsortedfirst(k_d, xj, 1, l_d, Base.Order.ForwardOrdering()) - 1
                end
            end
        end
        $ex
        dest
    end
end

function getindex{T,N,TCoefs,IT<:DimSpec{Gridded},K,P}(itp::GriddedInterpolation{T,N,TCoefs,IT,K,P}, x...)
    dest = Array(T, map(length, x))::Array{T,N}
    getindex!(dest, itp, x...)
end
