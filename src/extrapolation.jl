abstract ExtrapolationBehavior

immutable ExtrapError <: ExtrapolationBehavior end
function extrap_transform_x(::OnGrid, ::ExtrapError, N)
    quote
        @nexprs $N d->(1 <= real(x_d) <= size(itp,d) || throw(BoundsError()))
    end
end
function extrap_transform_x(::OnCell, ::ExtrapError, N)
    quote
        @nexprs $N d->(1//2 <= real(x_d) <= size(itp,d) + 1//2 || throw(BoundsError()))
    end
end

immutable ExtrapNaN <: ExtrapolationBehavior end
function extrap_transform_x(::OnGrid, ::ExtrapNaN, N)
    quote
        @nexprs $N d->(1 <= real(x_d) <= size(itp,d) || return convert(T, NaN))
    end
end
function extrap_transform_x(::OnCell, ::ExtrapNaN, N)
    quote
        @nexprs $N d->(1//2 <= real(x_d) <= size(itp,d) + 1//2 || return convert(T, NaN))
    end
end

immutable ExtrapConstant <: ExtrapolationBehavior end
function extrap_transform_x(::OnGrid, ::ExtrapConstant, N)
    quote
        @nexprs $N d->(x_d = clamp(x_d, 1, size(itp,d)))
    end
end
function extrap_transform_x(::OnCell, ::ExtrapConstant, N)
    quote
        @nexprs $N d->(x_d = clamp(x_d, 1//2, size(itp,d)+1//2))
    end
end

immutable ExtrapReflect <: ExtrapolationBehavior end
function extrap_transform_x(::OnGrid, ::ExtrapReflect, N)
    quote
        @nexprs $N d->begin
            # translate x_d to inside the domain, and count the translations
            ntransl = 0
            while x_d < 1
                x_d += size(itp, d) - 1
                ntransl += 1
            end
            while x_d > size(itp, d)
                x_d -= size(itp, d) - 1
                ntransl += 1
            end

            # if odd number of translations, also reflect inside the domain
            if ntransl > 0 && mod(ntransl, 2) != 0
                x_d = size(itp, d) + 1 - x_d
            end
        end
    end
end
function extrap_transform_x(::OnCell, ::ExtrapReflect, N)
    quote
        @nexprs $N d->begin
            # translate x_d to inside the domain, and count the translations
            ntransl = 0
            while x_d < .5
                x_d += size(itp,d)
                ntransl += 1
            end
            while x_d > size(itp,d) + .5
                x_d -= size(itp,d)
                ntransl += 1
            end
            
            # if odd number of translations, also reflect inside the domain
            if ntransl > 0 && mod(ntransl, 2) != 0
                x_d = size(itp, d) + 1 - x_d
            end 
        end
    end
end

immutable ExtrapPeriodic <: ExtrapolationBehavior end
function extrap_transform_x(::GridRepresentation, ::ExtrapPeriodic, N)
    :(@nexprs $N d->(x_d = mod1(x_d, size(itp,d))))
end

immutable ExtrapLinear <: ExtrapolationBehavior end
function extrap_transform_x(::OnGrid, ::ExtrapLinear, N)
    quote
        @nexprs $N d->begin
            if x_d < 1
                fx_d = x_d - convert(typeof(x_d), 1)

                k = -4*itp.coefs[1]
                return itp[1] - k * fx_d
            end
            if x_d > size(itp, d)
                s_d = size(itp,d)
                fx_d = x_d - convert(typeof(x_d), s_d)

                k = itp[s_d] - itp[s_d - 1]
                return itp[s_d] + k * fx_d
            end
        end
    end
end
extrap_transform_x(::OnCell, e::ExtrapLinear, N) = extrap_transform_x(OnGrid(), e, N)
