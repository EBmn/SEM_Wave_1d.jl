

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_2d_Updated.jl")
using .SEM_Wave_2d_Updated


include("MMS.jl")
using .MMS




function main()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######

    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0


    N = 8   # degree of interpolation polynomials
    Kx = 8 # number of elements in x-direction
    Ky = 8 # number of elements in y-direction

    heaviside(x) = 0.5 * (sign(x) + 1)

    c_square(x, y) = 1

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    Tend = 0.25
    nsteps = 200

    fVals = zeros(length(simul.y), length(simul.x))

    alphas = [0; 1/sqrt(2); 0; 1/sqrt(2)] #impedance up and down, Neumann left and right

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    omega = 7.1

    Tend = 5.0
    nsteps = 2500

    # construct initial data
    #uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    uStart = exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'
    #uStart = 2 .- (abs.(simul.y) .+ abs.(simul.x'))

    #uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    #uStartDer = 1e3*exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'

    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alpas, g, true, 10, 5)


end

