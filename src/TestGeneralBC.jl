

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_1d_new.jl")
using .SEM_Wave_1d


include("MMS.jl")
using .MMS



function main()

    ##########################
    #set up the time interval, the number of elements, the degree of the interpolation polynomials

    xl = 0.0
    xr = 1.0
    NumberOfNodes = 6
    #NumberOfNodes = 3211
    nodes = collect(LinRange(xl, xr, NumberOfNodes))

    N = 5 #number of points we interpolate in for each element
    
    simul = SEM_Wave_1d.SEM_Wave(nodes, N)

    simul.omega = 3.7
    # SEM_Wave_1d_new.SetSimulParams(Tend,nsteps,fVals,uData,uDataDer)
    simul.Tend = 10.0
    
    #delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    #simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
    simul.nsteps = 1000
    
    #println(simul.timestep)

    #set up some inital conditions and the values of the driving term f
    simul.fVals = 50*exp.(-50*(simul.x .- (xr-xl)/2).^2)

    #simul.fVals = zeros(length(simul.x))
    #simul.fVals[120:145] = 1500*ones(26)
    #simul.fVals = ones(length(simul.x))
    #simul.fVals = 3000*exp.(-((simul.x.-0.5)/0.05).^2)


    
    #simul.g = [1.0, 4.0]
    #simul.g = [1.0, -2.0]
    simul.g = [0.0, 0.0]


    #uStart = 5*exp.(-((simul.x.-0.5)/0.1).^2)
    
    #uStart = -ones(length(simul.x))
    uStart = zeros(length(simul.x))
    #uStart[5:7] = ones(3)
    #uStartDer = zeros(length(simul.x))
    #uStartDer = 80*exp.(-((simul.x.-0.5)/0.1).^2)
    uStartDer = zeros(length(simul.x))
    


    simul.bc = [0.0, 1.0]
    #delta = 0.28
    #simul.bc = [1/sqrt(2) + delta, sqrt(1/2 - sqrt(2)*delta - delta^2)]
    #simul.bc = [1/sqrt(2), 1/sqrt(2)]
    

    animate = false
    snapshotFrequency = 5

    #println(simul.M)
    #println(simul.M_b)
    
    SEM_Wave_1d.Simulate(simul, uStart, uStartDer, animate, snapshotFrequency)
    #SEM_Wave_1d.Simulate(simul, uStart, uStartDer, false)

    #plt = plot(simul.x, simul.uNow)

    println(simul.uNow)
    #savefig(plt, "solTest")

end

