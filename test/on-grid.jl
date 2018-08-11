module OnGridTests

using Interpolations, Test

nx, ny, nz = 10, 8, 9
xg, yg, zg = 1:nx, 1:ny, 1:nz

f1(x) = sin(2pi/nx * (x-1))
f2(x,y) = f1(x) * cos(4pi/ny * (y-1))
f3(x,y,z) = f2(x,y) * log(z)

simpledegs = (Linear,Constant)
bcdegs = (Quadratic,)
gridbehvs = (OnCell,OnGrid)
bcs = (Flat,Line,Free,Periodic,Reflect)

function OneD()
    println("Testing evaluation on grid and boundary in 1D...")
    for deg in simpledegs, gb in gridbehvs
        fg = Float64[f1(x) for x in xg]

        itp = Interpolation(fg, deg(gb()), ExtrapError())

        for i in xg
            @test ≈(itp[i],fg[i],atol=sqrt(eps()))
        end

        if gb == OnCell
            for xb in [.5, nx+.5]
                itp[xb] # if we're out of bounds, this throws
            end
        end
    end

    for deg in bcdegs, gb in gridbehvs, bc in bcs
        fg = Float64[f1(x) for x in xg]

        itp = Interpolation(fg, deg(bc(),gb()), ExtrapError())

        for i in xg
            @test ≈(itp[i],fg[i],atol=sqrt(eps()))
        end

        if gb == OnCell
            for xb in [.5, nx+.5]
                itp[xb]
            end
        end
    end
end

function TwoD()
    println("Testing evaluation on grid and boundary in 2D...")
    for deg in simpledegs, gb in gridbehvs
        fg = Float64[f2(x,y) for x in xg, y in yg]
        itp = Interpolation(fg, deg(gb()), ExtrapError())

        for i in xg, j in yg
            @test ≈(itp[i,j],fg[i,j],atol=sqrt(eps()))
        end

        if gb == OnCell
            for xb in [.5, nx+.5], yb in [.5,ny+.5]
                itp[xb,yb]
            end
        end
    end

    for deg in bcdegs, gb in gridbehvs, bc in bcs
        fg = Float64[f2(x,y) for x in xg, y in yg]
        itp = Interpolation(fg, deg(bc(),gb()), ExtrapError())

        for i in xg, j in yg
            @test ≈(itp[i,j],fg[i,j],atol=sqrt(eps()))
        end

        if gb == OnCell
            for xb in [.5, nx+.5], yb in [.5,ny+.5]
                itp[xb,yb]
            end
        end
    end
end

function ThreeD()
    println("Testing evaluation on grid and boundary in 3D...")
    for deg in simpledegs, gb in gridbehvs
        fg = Float64[f3(x,y,z) for x in xg, y in yg, z in zg]
        itp = Interpolation(fg, deg(gb()), ExtrapError())

        for i in xg, j in yg, k in zg
            @test ≈(itp[i,j,k],fg[i,j,k],atol=sqrt(eps()))
        end

        if gb == OnCell
            for xb in [.5, nx+.5], yb in [.5,ny+.5], zb in [.5, nz+.5]
                itp[xb,yb,zb]
            end
        end
    end

    for deg in bcdegs, gb in gridbehvs, bc in bcs
        fg = Float64[f3(x,y,z) for x in xg, y in yg, z in zg]
        itp = Interpolation(fg, deg(bc(),gb()), ExtrapError())

        for i in xg, j in yg, k in zg
            @test ≈(itp[i,j,k],fg[i,j,k],atol=sqrt(eps()))
        end

        if gb == OnCell
            for xb in [.5, nx+.5], yb in [.5,ny+.5], zb in [.5, nz+.5]
                itp[xb,yb,zb]
            end
        end
    end
end

OneD()
TwoD()
ThreeD()

end
