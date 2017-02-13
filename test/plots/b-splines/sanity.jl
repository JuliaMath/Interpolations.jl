using Interpolations
#using DataFrames, Gadfly

f(x) = sin(2pi * (x-1) / 10)
dfs = Array{DataFrame}(0)
xrange = .5:.1:11.5

for (i,it) in enumerate((
    BSpline(Constant()),
    BSpline(Linear()),
    BSpline(Quadratic(Free()))
    ))
    itp = interpolate(f(1:11), it, OnCell())
    push!(dfs, DataFrame(x=xrange,y=[itp[x]+(i-1)*3 for x in xrange],t="$it"))
end

push!(dfs, DataFrame(x=xrange, y=f(xrange),t="Data"))
df = vcat(dfs)
#display(plot(df,x=:x,y=:y,color=:t,Geom.path,Guide.colorkey("Dataset")))
    