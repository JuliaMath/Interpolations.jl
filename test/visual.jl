module VisualTests

print("loading packages...")
using Gadfly, Interpolations
println("done!")
const nx = 10
xg = 1:nx
xf = -.2nx:.1:4.3nx

f1(x) = sin(2pi/nx * (x-1))

y1 = map(f1,xg)

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
    )
    itp = Interpolation(y1, IT(), EB())

    display(plot(
        layer(x=xg,y=y1,Geom.point),
        layer(x=xf,y=[itp[x] for x in xf],Geom.path),
        Guide.title("$IT, $EB")
    ))
end

end
