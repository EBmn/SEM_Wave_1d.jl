

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_1d_new.jl")
using .SEM_Wave_1d


include("MMS.jl")
using .MMS




function main()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    xl = 0
    xr = 1

    delta = 0.2
    bc = [1/sqrt(2) + delta, sqrt(1/2 - sqrt(2)*delta - delta^2)]

    #bc = [0.0, 1.0]
    #bc = [1/sqrt(2), 1/sqrt(2)]
    
    
    numberOfNodes = 24


    N = 4 # number of points we interpolate in for each element

    plt = plot()

    nodes = collect(LinRange(xl, xr, numberOfNodes))

    simul = SEM_Wave_1d.SEM_Wave(nodes, N)

    ###### set up the problem ######
    #fVals = exp.(-50*(simul.x .- (xr-xl)/2).^2)
    #fVals = exp.(0*(simul.x .- (xr-xl)/2).^2)
    

    

    fVals = cos.(3*pi*simul.x) + 5*exp.(-50*(simul.x .- (xr-xl)/2).^2) + 7*exp.(-70*(simul.x .- (xr-xl)/1.2).^2)

    simul.g = [0.0, 0.0]
    omega = 1.27
    #simul.omega = 5*pi^2
    simul.Tend = 2*pi/omega
    #simul.fVals = (9*pi^2 - simul.omega^2)*cos.(3*pi*simul.x)
    #simul.fVals = cos.(3*pi*simul.x)
    
    
    #println(simul.nsteps)
    #println("CFL number: " * string(simul.timestep/delta_x))
    
    simul.bc = bc

    ###### compute an approximate Helmholtz solution ######
    tol = 1e-9
    sol1, nIter1 = SEM_Wave_1d.Waveholtz(simul, omega, fVals, tol)


    sol2, nIter2 = SEM_Wave_1d.WaveholtzGMRES(simul, omega, fVals, tol)
    
    println()
    println(nIter1)
    println(nIter2)
    println()
    #SEM_Wave_1d.Waveholtz(simul, 100)
    #SEM_Wave_1d.WaveholtzAnimation(simul, omega, fVals, 100)

    
    # plot the result
    #plt = plot!(simul.x, sol1, label=:false)
    #plot!(simul.x, sol2[1:length(simul.x)], label=:false)
    #plt = plot!(simul.x, simul.uFiltered -solution[1:length(simul.x)], label=:false)

    
    
    #code for verifying that the solution is correct:
    
    sol1_xx = SEM_Wave_1d.LaplaceTerm(simul, sol1) / simul.timestep^2

    sol2 = sol2[1:length(simul.x)]
    sol2_xx = SEM_Wave_1d.LaplaceTerm(simul, sol2) / simul.timestep^2
    
    

    plt = plot!(simul.x[2:end-1], sol1_xx[2:end-1] + omega^2*sol1[2:end-1] - fVals[2:end-1], label="fAppx1 - fVals")
    plot!(simul.x[2:end-1], sol2_xx[2:end-1] + omega^2*sol2[2:end-1] - fVals[2:end-1], label="fAppx2 - fVals")
    

    #plot!(simul.x[2:end-1], sol1_xx[2:end-1], label="2nd der")
    #plot!(simul.x[2:end-1], sol2_xx[2:end-1], label="2nd der")
    #plot!(simul.x[2:end-1],omega^2*sol1[2:end-1], label="u")
    
    savefig(plt, "WaveholtzTest")

end

