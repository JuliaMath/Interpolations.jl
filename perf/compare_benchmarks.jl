# NOTE: create symlinks called "results_old.json" and "results_new.json" in this directory
# to the results files you want to compare

using BenchmarkTools, PyPlot, Unitful

const tref = 1000

function xycmp(results_old, results_new, filterstrs...)
    x, y = Float64[], Float64[]
    for (k, v) in results_new
        passes = true
        for str in filterstrs
            passes &= occursin(str, k)
        end
        passes || continue
        if haskey(results_old, k)
            push!(x, minimum(results_old[k]).time/tref)
            push!(y, minimum(v).time/tref)
        end
    end
    x, y
end

function plotcmp(results_old, results_new, keys, filterstrs...)
    for key in keys
        results_old = results_old[key]
        results_new = results_new[key]
    end
    x, y = xycmp(results_old, results_new, filterstrs...)
    maxxy = max(maximum(x), maximum(y))
    scatter(x, y)
    plot([0, maxxy], [0, maxxy], "--")
    # title(string(keys, filterstrs))
    title(string(keys))
end

function plotcmppanels(results_old, results_new, keys, filterstrs...)
    nrows = ceil(Int, sqrt(length(keys)))
    ncols = ceil(Int, length(keys)/nrows)
    figure()
    for k = 1:length(keys)
        subplot(nrows, ncols, k)
        plotcmp(results_old, results_new, keys[k], filterstrs...)
    end
    suptitle(string(filterstrs))
end

results_old = BenchmarkTools.load("results_old.json")[1]
results_new = BenchmarkTools.load("results_new.json")[1]

results_old = results_old["bsplines"]
results_new = results_new["bsplines"]

for filt in (("use",), ("construct",))
    # for keys in (("constant",), ("linear",), ("quadratic",), ("cubic",))
    #     plotcmp(results_old, results_new, keys, filt...)
    # end
    keys = (("constant",), ("linear",), ("quadratic",), ("cubic",))
    plotcmppanels(results_old, results_new, keys, filt...)
end
