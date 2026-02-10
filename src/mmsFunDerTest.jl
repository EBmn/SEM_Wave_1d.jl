

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
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0
  

    

    N = 3   # degree of interpolation polynomials
    Kx = 4 # number of elements in x-direction
    Ky = 4 # number of elements in y-direction

    c_square(x, y) = 1 

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    Tend = 0.5
    nsteps = 200

    simul.useMMS = true
    
    simul.MMS_j.type = [3 3 3]

    simul.MMS_j.coeff[1, 1] = 1
    simul.MMS_j.coeff[1, 2] = 0.0
    simul.MMS_j.coeff[1, 3] = 5.0
    
    simul.MMS_j.coeff[2, 1] = 1
    simul.MMS_j.coeff[2, 2] = 0.0
    simul.MMS_j.coeff[2, 3] = 0.0

    simul.MMS_j.coeff[3, 1] = 1
    simul.MMS_j.coeff[3, 2] = 1
    simul.MMS_j.coeff[3, 3] = 0

    simul.MMS_c.type = [3 3]

    simul.MMS_c.coeff[1, 1] = 1
    simul.MMS_c.coeff[1, 2] = 0
    simul.MMS_c.coeff[1, 3] = 0

    simul.MMS_c.coeff[2, 1] = 1
    simul.MMS_c.coeff[2, 2] = 0.0
    simul.MMS_c.coeff[2, 3] = 1


    fVals = zeros(length(simul.y), length(simul.x))

    a = 0.23
    b = sqrt(1-a^2)
    bc = [a, b]

    g = zeros(length(simul.y), length(simul.x))


    simul.g = g
    alpha = simul.bc[1]
    beta = simul.bc[2]

    u_xx = zeros(length(simul.y), length(simul.x))
    u_yy = zeros(length(simul.y), length(simul.x))
    comparison = zeros(length(simul.y), length(simul.x))


    println("u_xx: " * string(MMS.MMSfun(-0.638, 0.2, 0.0, 2, 0, 0, simul.MMS_j)) * ", u_yy: " * string(MMS.MMSfun(0.7, 0.362, 0.0, 0, 2, 0, simul.MMS_j)))
    #println(MMS.MMSfun(0.3, 0.362, 0.0, 0, 0, 0, simul.MMS_j))


    for i = 1:length(simul.y)
        for j = 1:length(simul.x)
            #u_xx[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 2, 0, 0, simul.MMS_j)
            #u_yy[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 2, 0, simul.MMS_j)
            #comparison[i, j] = 20*simul.y[end+1-i].^3
        end
    end

    #plt = surface(simul.x, simul.y[end:-1:1], u_yy-comparison)
    #savefig(plt, "u_yyTester")

end

