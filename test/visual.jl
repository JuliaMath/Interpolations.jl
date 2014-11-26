module VisualTests

using Gadfly, Interpolations

const nx = 10
xg = 1:nx
xf = -2nx:.1:3nx

f1(x) = sin(2pi/nx * (x-1))

y1 = map(f1,xg)

for EB in (
        ExtrapNaN,
        ExtrapConstant,
        ExtrapReflect
    )

    for IT in (
            Constant{OnCell},
            Constant{OnGrid},
            Linear{OnCell},
            Linear{OnGrid},
            Quadratic{ExtendInner,OnCell},
            Quadratic{ExtendInner,OnGrid},
            Quadratic{Flat,OnCell},
            Quadratic{Flat,OnGrid},
        )
    
        itp = Interpolation(y1, IT(), EB())

        display(plot(
            layer(x=xg,y=y1,Geom.point),
            layer(x=xf,y=[itp[x] for x in xf],Geom.path),
            Guide.title("$(typeof(Interpolations.degree(IT()))), $(typeof(Interpolations.gridrepresentation(IT()))), $EB")
        ))
    end
end

for IT in (
        Constant{OnCell},
        Constant{OnGrid},
        Linear{OnCell},
        Linear{OnGrid},
        Quadratic{ExtendInner,OnCell},
        Quadratic{ExtendInner,OnGrid},
        Quadratic{Flat,OnCell},
        Quadratic{Flat,OnGrid},
    )

    # Treat ExtrapError specially, since it will throw bounds
    # errors for any x outside 1:nx
    itperr = Interpolation(y1, IT(), ExtrapError())
    if Interpolations.gridrepresentation(IT()) == OnCell()
        display(plot(
            layer(x=xg,y=y1,Geom.point),
            layer(x=.5:.1:nx+.5,y=[itperr[x] for x in .5:.1:nx+.5],Geom.path),
            Guide.title("$(Interpolations.degree(IT())), $(Interpolations.gridrepresentation(IT())), ExtrapError")
        ))
    else
        display(plot(
            layer(x=xg,y=y1,Geom.point),
            layer(x=1:.1:nx,y=[itperr[x] for x in 1:.1:nx],Geom.path),
            Guide.title("$(Interpolations.degree(IT())), $(Interpolations.gridrepresentation(IT())), ExtrapError")
        ))
    end
end

end