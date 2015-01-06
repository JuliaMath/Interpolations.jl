module VisualTests

using Gadfly, Interpolations
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
        (Quadratic{Line,OnCell}, ExtrapLinear),
        (Quadratic{Line,OnGrid}, ExtrapLinear),
        (Quadratic{Line,OnGrid}, ExtrapReflect),
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
    title = "$(IT.name.name){$(join(map(t -> t.name.name, IT.parameters), ','))}, $(EB.name.name)"
    push!(stuff, Guide.title(title))
    display(plot(stuff...))
end

nx = 10
xg = 1:nx
ny = 6
yg = 1:ny

f2(x,y) = sin(2pi/nx * (x-1)) * cos(2pi/ny * (y-1))
zg = Float64[f2(x,y) for x in xg, y in yg]

xf = -1:.1:nx+1
yf = -1:.1:ny+1

itp2 = Interpolation(zg, Quadratic(Flat(),OnCell()),ExtrapLinear())

display(plot(
    layer(x=xf,y=yf,z=[itp2[x,y] for x in xf, y in yf], Geom.contour),
    Guide.title("Quadratic{Flat,Oncell}, ExtrapLinear")
))

end

end