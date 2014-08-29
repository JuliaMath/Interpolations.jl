using Base.Test
import Interpolations

if !(isdefined(Main, :Gadfly) && isdefined(Main, :DataFrames))
	println("For a visualization of the tests, do")
	println(" using Gadfly, DataFrames")
	println(" include(joinpath(Pkg.dir(\"Interpolations\"), \"src/Interpolations.jl\"))")
	println(" include(joinpath(Pkg.dir(\"Interpolations\"), \"test\/runtests.jl\"))")
	println("preferrably in an IJulia notebook (it produces a lot of plots...)")
end
include("linear.jl")
