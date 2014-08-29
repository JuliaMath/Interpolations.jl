## Tests for linear interpolation with various extrapolation behaviors

## Construct ground-truth values for 1d interpolation
## To verify these data sets, do something like this:
# using Gadfly, DataFrames, Base.Test, Interpolations
# include(joinpath(Pkg.dir("Interpolations"), "test/linear.jl"))

f(x) = sin((x-3)*2pi/9 - 1)
xmax = 10

function A1linstep(step)
    ixs = step .<= xpos .< step+1
    dydx = f(step+1)-f(step)
    Astep = f(step) + (xpos[ixs]-step)*dydx
end

xpos = [-1.5xmax:0.1:xmax+2.5]
inbounds = 1 .<= xpos .<= xmax
A1inbounds = vcat([A1linstep(s) for s in 1:xmax-1]..., f(xmax))

A1 = similar(xpos)
A1[inbounds] = A1inbounds

A1error = fill(Inf, length(xpos))
A1error[inbounds] = A1inbounds  # use Inf as a sentinel for BoundsError

A1nan = fill(NaN, length(xpos))
A1nan[inbounds] = A1inbounds

# A1specified = copy(A1)
# A1specified[xpos .< 1] = A1specified[xpos .> xmax] = -1

A1constant = copy(A1)
A1constant[xpos .< 1] = A1inbounds[1]
A1constant[xpos .> xmax] = A1inbounds[end]

l1 = last(find(xpos .< 1)); 
l2 = findfirst(xpos .> xmax)-1

A1linear = copy(A1)
dydx1 = f(2)-f(1)
A1linear[1:l1] = f(1) + dydx1 * (xpos[1:l1]-1)
dydx2 = f(xmax)-f(xmax-1)
A1linear[l2:end] = f(xmax) + dydx2 * (xpos[l2:end]-xmax)


# Plot the values, to make sure we've done it right
if isdefined(Main, :Gadfly) && isdefined(Main, :DataFrame)
    Aall = vcat(
        # DataFrame(x = xpos, y = A1na,  eb = "BCna"),
        DataFrame(x = xpos[inbounds], y = A1error[inbounds], eb = "Function values"), # plot only inbounds values

        DataFrame(x = xpos[!inbounds], y = A1constant[!inbounds], eb = "ExtrapConstant"),
        DataFrame(x = xpos[!inbounds], y = A1linear[!inbounds], eb = "ExtrapLinear"),
        # DataFrame(x = xpos[!inbounds], y = A1reflect[!inbounds], eb = "ExtrapReflect"),
        # DataFrame(x = xpos[!inbounds], y = A1reflect2[!inbounds], eb = "ExtrapReflect2"),
        DataFrame(x = xpos[!inbounds], y = A1periodic[!inbounds], eb = "ExtrapPeriodic"),
        #DataFrame(x = xpos[!inbounds], y = A1specified[!inbounds], eb = "ExtrapSpecified"),
    )
    set_default_plot_size(20cm, 10cm)
    display(plot(Aall,
        layer(x = :x,
        y = :y,
        color = :eb,
        Geom.point, 
        #Guide.xticks(ticks=1:4)
        ),
        layer(x=xpos[inbounds],y=f(xpos[inbounds]),Geom.path),
        Guide.colorkey("Extrap behavior"),
        Scale.y_continuous(minvalue=-2,maxvalue=2)
    ))
end

# for (EB,correct) in ((Interpolations.ExtrapError,A1error), (Interpolations.ExtrapNaN,A1nan), (Interpolations.ExtrapConstant, A1constant))
#     G = Interpolations.Interpolation(f(1:xmax), Interpolations.Linear, EB)
#     if isdefined(Main, :Gadfly) 
#         if EB == Interpolations.ExtrapPeriodic
#             display(plot(
#                 layer(x=xpos[inbounds],y=f(xpos[inbounds]), Geom.path),
#                 layer(x=xpos,y=correct,Geom.point,Theme(default_color=color("green"))),
#                 layer(x=xpos,y=[G[x] for x in xpos],Geom.point, Theme(default_color=color("red"))),
#                 Guide.xticks(ticks=[-1:xmax+2]),
#                 Guide.title("$EB")
#             ))
#         end
#     else
#         println(EB)
#     end
#     for i = 1:length(xpos)
#         x = xpos[i]
#         y = correct[i]
#         if !isinf(y)
#             try
#                 @test_approx_eq y G[x]
#             catch err
#                 @show x, y, G[x]
#                 rethrow(err)
#             end
#         else
#             @test_throws BoundsError G[x]
#         end
#     end
# end
