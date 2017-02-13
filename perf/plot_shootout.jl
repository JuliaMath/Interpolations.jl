using Gadfly, DataFrames

i = indexin([Float64], [eltypes...])[1]
tconstr = tconstr[:,:,i]
teval = teval[:,:,i]

# Convert to dataframes, and for teval use points/second
dim = ones(Int,size(tconstr,1))*(1:size(tconstr,2))'
namelist = collect(take(cycle(testnames), length(dim)))
ordlist =  collect(take(cycle(testord), length(dim)))
dfconstr = DataFrame(t = vec(tconstr), name = namelist, ord = ordlist, dim = vec(dim))
dfeval   = DataFrame(rate = vec(10^6./teval), name = namelist, ord = ordlist, dim = vec(dim))

flagBspline = (dfconstr[:name] .== "IBSpline") | (dfconstr[:name] .== "Grid")
dfBconstr = dfconstr[flagBspline, :]
dfBeval   = dfeval[flagBspline, :]
dfGconstr = dfconstr[!flagBspline, :]
dfGeval   = dfeval[!flagBspline, :]

yscale = Scale.y_continuous

# Plot the B-spline methods
plotsize = (6inch, 8inch)
set_default_plot_size(plotsize...)
pBc = plot(dfBconstr, xgroup="dim", ygroup="ord", x="name", y="t", Geom.subplot_grid(Geom.point), yscale, Guide.ylabel("Time (s) by order"))
draw(PNG("constructionB.png", plotsize...), pBc)

pBr = plot(dfBeval, xgroup="dim", ygroup="ord", x="name", y="rate", Geom.subplot_grid(Geom.point), yscale, Guide.ylabel("Throughput (pts/s) by order"))
draw(PNG("rateB.png", plotsize...), pBr)

# Plot the gridded methods
plotsize = (10inch, 8inch)
set_default_plot_size(plotsize...)

pGc = plot(dfGconstr, xgroup="dim", ygroup="ord", x="name", y="t", Geom.subplot_grid(Geom.point), yscale, Guide.ylabel("Time (s) by order"))
draw(PNG("constructionG.png", plotsize...), pGc)

pGr = plot(dfGeval, xgroup="dim", ygroup="ord", x="name", y="rate", Geom.subplot_grid(Geom.point), Scale.y_sqrt, Guide.ylabel("Throughput (pts/s) by order"))
draw(PNG("rateG.png", plotsize...), pGr)