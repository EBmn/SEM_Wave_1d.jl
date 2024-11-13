
using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using TimerOutputs


include("SEM_Wave_1d.jl")

#=
    Produces convergence plots SEM_Wave for a fixed polynomial degree and increasing number of elements
    We let the final run be the reference solution, and compare more coarse runs to it.
    Finally, we plot these differences in a loglog plot, hopefully finding reasonable convergence rates
=#

function main()

    numberOfRuns = 8
    numberOfNodesList = zeros(numberOfRuns)
    timestepsList = zeros(numberOfRuns)

    for j = 1:numberOfRuns
        numberOfNodesList[j] = 2^j + 1
        timestepsList[j] = 20*2^j + 1
    end

    numberOfNodesList = Int.(numberOfNodesList)
    timestepsList = Int.(timestepsList)

    println(timestepsList)

    ordersList = [3, 5, 6, 7]
    
    omega = 6.0

    GenerateConvPlot(ordersList, numberOfNodesList, timestepsList, numberOfRuns, omega)
    
end


function GenerateConvPlot(ordersList, numberOfNodesList, timestepsList, numberOfRuns, omega)

    plt = plot()

    for n = 1:length(ordersList)
        #the number of quadrature points per element; frequency parameter
        N = ordersList[n]
        
        #number of timesteps + outer limit of time interval
        Tend = 1.0
        #nsteps = timestepsList[n]

        #outer limits of the interval we are solving the equation on
        xl = 0
        xr = 1

        #where we store the difference between the jth run and the final one
        diffs = zeros(numberOfRuns-1)
        reference = zeros((N-1)*(numberOfNodesList[end]-1) + 1)

        for j = numberOfRuns:-1:1
            #solve the wave equation with the different numbers of nodes
            
            #set up some inital conditions and the values of the driving term f

            nsteps = timestepsList[j]

            numberOfNodes = numberOfNodesList[j]
            
            nodes = zeros(numberOfNodes)
            
            range = LinRange(xl, xr, numberOfNodes)
            
            for i = 1:numberOfNodes
                nodes[i] = range[i]
            end

            quadpoints, ~ = gausslobatto(N)
            totalPointCount = (length(nodes)-1) * (N-1) + 1

            x = SEM_Wave_1d.ConstructX(nodes, quadpoints)
            fVals = zeros(totalPointCount)
            uStart = zeros(totalPointCount)
            uStartDer = zeros(totalPointCount)

            
            #the driving term is a function that is C^{\infty}_0 on our interval, supported in {x : |centreOfSupport - x | \leq 2*epsilon]}
            centreOfSupport = (xr - xl)/2
            epsilon = 0.1

            #fill fVals with appropriate values
            for i = 1:totalPointCount
                fVals[i] = 50*exp(-1/(1-min(1, (x[i] - centreOfSupport)^2/(2*epsilon)^2)))
            end

            animate = false
            snapshotFrequency = 0   #as animate = false, this value is not used anywhere

            output = SEM_Wave_1d.Simulate(nodes, nsteps, Tend, N, fVals, omega, uStart, uStartDer, animate, snapshotFrequency)

            if (j == numberOfRuns)
                reference = output
            else
                diffs[j] = DiffCalc(output, reference, j, numberOfRuns)
            end

        end

        println(diffs)

        plot!(numberOfNodesList[1:end-1], diffs, xscale=:log10, yscale=:log10, xlims=(numberOfNodesList[1], numberOfNodesList[end]), label = "N = " * string(N), legend=:bottomleft)            
        println("done with N = " * string(N))        

    end

    display(plt)

end


function DiffCalc(output, reference, j, numberOfRuns)
    #returns the maximum of the difference between the values of output and reference in the nodes that they share
    sharedNodeIndex = 0
    sharedNodes = zeros(length(output))

    for i = 1:length(sharedNodes)
        sharedNodeIndex = 1 + (i-1)*2^(numberOfRuns-j)
        sharedNodes[i] = reference[sharedNodeIndex]
    end

    return maximum(abs.(sharedNodes - output))

end

