


using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_2d.jl")
using .SEM_Wave_2d


include("MMS.jl")
using .MMS




function main()


    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######


    # test 1
    #=
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0
    =#

    # test 2
    xl = -1.0
    xr = 2.0
    yl = -2.0
    yr = 1.0


    N = 4   # degree of interpolation polynomials
    Kx = 16 # number of elements in x-direction
    Ky = 16 # number of elements in y-direction

    # test 1
    #c_square(x, y) = 1.0 

    # test 2
    c_square(x, y) = 1.0 + x^2 + 2*y^2
    


    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    
    # test 1
    #=
    p = 1
    q = 4
    s = 3
    t = 2
    =#

    # test 2
    p = 1
    q = 0
    s = 1
    t = 0

    U = (simul.y).^q * (simul.x').^p
    V = (simul.y).^t * (simul.x').^s

    # u = x^py^q and v = x^sy^t gives a  gradTntegral over [0, 1]^2 of ps/((p+s-1)(q+t+1)) + qt/((p+s+1)(q+t-1))

    gradIntVal = SEM_Wave_2d.GradIntegral(simul, U, V)
    
    # test 1
    #trueGradIntVal = p*s/((p+s-1)*(q+t+1)) + q*t/((p+s+1)*(q+t-1)) 

    # test 2
    trueGradIntVal = 36



    println("should be " * string(trueGradIntVal) * ": " * string(gradIntVal))
    println("error: " * string(abs(trueGradIntVal - gradIntVal)))

end

