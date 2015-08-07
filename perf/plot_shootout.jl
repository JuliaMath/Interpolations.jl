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

set_default_plot_size(15cm, 15cm)

pc = plot(dfconstr, xgroup="dim", ygroup="ord", x="name", y="t", Geom.subplot_grid(Geom.point), Scale.y_sqrt, Guide.ylabel("Time (s) by order"))
draw(PNG("construction.png", 8inch, 6inch), pc)

pr = plot(dfeval, xgroup="dim", ygroup="ord", x="name", y="rate", Geom.subplot_grid(Geom.point), Scale.y_sqrt, Guide.ylabel("Throughput (pts/s) by order"))
draw(PNG("rate.png", 8inch, 6inch), pr)
