

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_1d_cleaned.jl")
using .SEM_Wave_1d


include("MMS.jl")
using .MMS




function main()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    xl = 0.0
    xr = 5.0

    delta = 0.25
    bc = [1/sqrt(2) + delta, sqrt(1/2 - sqrt(2)*delta - delta^2)]

    #bc = [0.0, 1.0]
    bc = [1/sqrt(2), 1/sqrt(2)]
    

    N = 1 # degree of interpolation polynomials
    K = 150 # number of elements

    #c_square(x) = 1
    #c_square(x) = 1 - x.*(x-5)./2
    c_square(x) = 2 + sign(x-4)

    simul = SEM_Wave_1d.SEM_Wave([xl, xr], K, N, c_square)
    

    plt = plot()

    ###### set up the problem ######
    fVals = 50*exp.(-50*(simul.x .- (xr-xl)/2).^2)
    #fVals = cos.(3*pi*simul.x) + 5*exp.(-50*(simul.x .- (xr-xl)/2).^2) + 7*exp.(-70*(simul.x .- (xr-xl)/1.2).^2)
    #fVals = zeros(length(simul.x))
    

    

    g = [0.0, 0.0]
    omega = 28.7
    
    Tend = 10.0
    nsteps = 800

    uStart = 8*exp.(-((simul.x.-0.5)/0.1).^2)
    uStartDer = zeros(length(simul.x))
    
    
    #println(simul.nsteps)
    #println("CFL number: " * string(simul.timestep/delta_x))


    #=

    ###### compute an approximate Helmholtz solution ######
    tol = 1e-12
    sol1, nIter1 = SEM_Wave_1d.Waveholtz(simul, fVals, omega, bc, g, tol)


    sol2, nIter2 = SEM_Wave_1d.WaveholtzGMRES(simul, fVals, omega, bc, g, tol)

    savefig(plt, "WaveholtzTest")
    
    #println()
    #println(nIter1)
    #println(nIter2)
    #println()
    #SEM_Wave_1d.Waveholtz(simul, 100)
    #SEM_Wave_1d.WaveholtzAnimation(simul, omega, fVals, 100)

    
    # plot the result
    plt = plot!(simul.x, sol1, label=:false)
    plot!(simul.x, sol2[1:length(simul.x)], label=:false)
    #plt = plot!(simul.x, simul.uFiltered -solution[1:length(simul.x)], label=:false)

    
    #=
    #code for verifying that the solution is accurate:
    
    sol1_xx = SEM_Wave_1d.LaplaceTerm(simul, sol1) / simul.timestep^2

    sol2 = sol2[1:length(simul.x)]
    sol2_xx = SEM_Wave_1d.LaplaceTerm(simul, sol2) / simul.timestep^2
    
    

    

    plt = plot!(simul.x[2:end-1], sol1_xx[2:end-1] + omega^2*sol1[2:end-1] - fVals[2:end-1], label="fAppx1 - fVals")
    plot!(simul.x[2:end-1], sol2_xx[2:end-1] + omega^2*sol2[2:end-1] - fVals[2:end-1], label="fAppx2 - fVals")
    

    #plot!(simul.x[2:end-1], sol1_xx[2:end-1], label="2nd der")
    #plot!(simul.x[2:end-1], sol2_xx[2:end-1], label="2nd der")
    #plot!(simul.x[2:end-1],o mega^2*sol1[2:end-1], label="u")
    
    savefig(plt, "WaveholtzTest")
    =#    
    
    =#
    SEM_Wave_1d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, true, 2)
    #println(simul.uNow)


    

end

