module VisualTests

print("loading packages...")
using Gadfly, Interpolations
println("done!")
const nx = 10
xg = 1:nx
xf = -1.2nx:.1:3.3nx

f1(x) = sin(2pi/nx * (x-1))

y1 = map(f1,xg)

p = plot()

if true

for (IT, EB) in (
        (Constant{OnCell}, ExtrapConstant),
        (Constant{OnGrid}, ExtrapConstant),
        (Linear{OnCell}, ExtrapNaN),
        (Linear{OnGrid}, ExtrapNaN),
        (Quadratic{Flat,OnCell}, ExtrapConstant),
        (Quadratic{Flat,OnGrid}, ExtrapConstant),
        (Quadratic{LinearBC,OnCell}, ExtrapLinear),
        (Quadratic{LinearBC,OnGrid}, ExtrapLinear),
        (Quadratic{LinearBC,OnGrid}, ExtrapReflect),
        (Quadratic{Flat,OnCell}, ExtrapReflect),
        (Quadratic{Free,OnCell}, ExtrapReflect),
        (Quadratic{Free,OnGrid}, ExtrapReflect),
        (Quadratic{Periodic,OnGrid}, ExtrapPeriodic),
        (Quadratic{Periodic,OnCell}, ExtrapPeriodic),
    )
    itp = Interpolation(y1, IT(), EB())

    stuff = Any[]

    push!(stuff, layer(x=xg,y=y1,Geom.point,Theme(default_color=color("green"))))
    push!(stuff, layer(x=xf,y=[itp[x] for x in xf],Geom.path,Theme(default_point_size=2px)))
    push!(stuff, Guide.title("$IT, $EB"))
    if boundarycondition(IT()) == Periodic()
        push!(stuff, layer(x=xg,y=y1,Geom.point,Theme(default_color=color("green"))))
    end
    display(plot(stuff...))
end

end

end