module VisualTests

using Gadfly, Interpolations
const nx = 10
xg = 1:nx
xf = -1.2nx:.1:3.3nx

f1(x) = sin(2pi/nx * (x-1))

y1 = f1.(xg)

p = plot()

if true

Btypes = (Periodic, Flat, Line, Free, Reflect)
Itypes = (
    Constant, Linear,
    [Quadratic{T} for T in Btypes]...,
    [Cubic{T} for T in Btypes[1:end-1]]...,  # no Reflect for Cubic
)
Etypes = (Flat, Linear, Reflect, Periodic)
Gtypes = (OnCell, OnGrid)

for IT in Itypes, GT in Gtypes, ET in Etypes
    itp = extrapolate(interpolate(y1, BSpline(IT()), GT()), ET())

    stuff = Any[]

    push!(stuff, layer(x=xg,y=y1,Geom.point,Theme(default_color=colorant"green")))
    push!(stuff, layer(x=xf,y=[itp[x] for x in xf],Geom.path,Theme(point_size=2px)))
    title = "$(IT.name.name){$(join(map(t -> t.name.name, IT.parameters), ','))}, $(ET.name.name)"
    push!(stuff, Guide.title(title))
    display(plot(stuff...))
end

nx = 10
xg = 1:nx
ny = 6
yg = 1:ny

f2(x,y) = sin(2pi/nx * (x-1)) * cos(2pi/ny * (y-1))
zg = f2.(xg, yg')

xf = -1:.1:nx+1
yf = -1:.1:ny+1

itp2 = extrapolate(interpolate(zg, BSpline(Quadratic{Flat}()), OnCell()), Linear())

display(plot(
    layer(x=xf,y=yf,z=[itp2[x,y] for x in xf, y in yf], Geom.contour),
    Guide.title("Quadratic{Flat}, Oncell, Linear")
))

end

end
