
module OnGridTests

using Interpolations, Base.Test

nx, ny, nz = 10, 8, 9
xg, yg, zg = 1:nx, 1:ny, 1:nz

f1(x) = sin(2pi/nx * (x-1))
f2(x,y) = f1(x) * cos(4pi/ny * (y-1))
f3(x,y,z) = f2(x,y) * log(z)

simpledegs = (Linear,Constant)
bcdegs = (Quadratic,)
gridbehvs = (OnCell,OnGrid)
bcs = (Flat,Line,Free,Periodic)

function OneD()
    println("Testing on-grid evaluation in 1D...")
    for deg in simpledegs, gb in gridbehvs
        fg = Float64[f1(x) for x in xg]
        
        itp = Interpolation(fg, deg(gb()), ExtrapError())

        for i in xg
            @test_approx_eq_eps itp[i] fg[i] sqrt(eps())
        end
    end

    for deg in bcdegs, gb in gridbehvs, bc in bcs
        fg = Float64[f1(x) for x in xg]

        itp = Interpolation(fg, deg(bc(),gb()), ExtrapError())

        for i in xg
            @test_approx_eq_eps itp[i] fg[i] sqrt(eps())
        end
    end
end

function TwoD()
    println("Testing on-grid evaluation in 2D...")
    for deg in simpledegs, gb in gridbehvs
        fg = Float64[f2(x,y) for x in xg, y in yg]
        itp = Interpolation(fg, deg(gb()), ExtrapError())

        for i in xg, j in yg
            @test_approx_eq_eps itp[i,j] fg[i,j] sqrt(eps())
        end
    end

    for deg in bcdegs, gb in gridbehvs, bc in bcs
        fg = Float64[f2(x,y) for x in xg, y in yg]
        itp = Interpolation(fg, deg(bc(),gb()), ExtrapError())

        for i in xg, j in yg
            @test_approx_eq_eps itp[i,j] fg[i,j] sqrt(eps())
        end
    end
end

function ThreeD()
    println("Testing on-grid evaluation in 3D...")
    for deg in simpledegs, gb in gridbehvs
        fg = Float64[f3(x,y,z) for x in xg, y in yg, z in zg]
        itp = Interpolation(fg, deg(gb()), ExtrapError())

        for i in xg, j in yg, k in zg
            @test_approx_eq_eps itp[i,j,k] fg[i,j,k] sqrt(eps())
        end
    end

    for deg in bcdegs, gb in gridbehvs, bc in bcs
        fg = Float64[f3(x,y,z) for x in xg, y in yg, z in zg]
        itp = Interpolation(fg, deg(bc(),gb()), ExtrapError())

        for i in xg, j in yg, k in zg
            @test_approx_eq_eps itp[i,j,k] fg[i,j,k] sqrt(eps())
        end
    end
end

OneD()
TwoD()
ThreeD()

end
