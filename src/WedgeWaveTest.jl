


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
    

    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0



    N = 5   # degree of interpolation polynomials
    Kx = 150 # number of elements in x-direction
    Ky = 150 # number of elements in y-direction

    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    omega = 4*pi
    a = 1/sqrt(2)
    b = sqrt(1-a^2)
    bc = [a, b]
    Tend = 2.5

    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))

    nsteps = Integer(ceil(1.5*Tend * cMax * (1/delta_x + 1/delta_y)))
    nsteps = 2*nsteps

    println("timestep is " * string(Tend/nsteps))

    delta(x) = Float64(x==0.0) # the best implementation known to man
    f(x, y) = omega^2 * delta(abs(x - 300.0))' * delta(abs(y))
    #fVals = f.(simul.x', simul.y)
    fVals = omega^2 * exp.(-0.01*(simul.y .- 20.0).^2) * exp.(-0.01*(simul.x .- 300.0).^2)'

    g = zeros(length(simul.y), length(simul.x))


    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))

    useMMS = false
    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, true, 1, 2)
    println(maximum(simul.uNow))
    

end

