
using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using TimerOutputs


include("SEM_Wave_1d.jl")

#=
Here we conduct an MMS-test for the Laplacian and plot the results for different values of N
=#
function main()

    numberOfRuns = 16
    numberOfNodesList = zeros(numberOfRuns)

    for j = 1:numberOfRuns
        #numberOfNodesList[j] = 5*2^j + 1
        numberOfNodesList[j] = 2^j + 1
    end

    numberOfNodesList = Int.(numberOfNodesList)

    ordersList = [2, 3, 4, 5, 6, 7]

    GenerateConvPlot(ordersList, numberOfNodesList, numberOfRuns)

end



function GenerateConvPlot(ordersList, numberOfNodesList, numberOfRuns)

    plt = plot()

    for n = 1:length(ordersList)
        #the number of quadrature points per element; frequency parameter
        N = ordersList[n]

        #just to have something to initialise a SEM_Wave-struct with
        nsteps = 1
        Tend = 1.0
        omega = 1.0
        
        #outer limits of the interval we are considering
        xl = 0
        xr = 1

        #where we store the difference between the jth run and the value found symbolically
        diffs = zeros(numberOfRuns)
        reference = zeros((N-1)*(numberOfNodesList[end]-1) + 1)

        for j = 1:numberOfRuns
            #solve the wave equation with the different numbers of nodes
            
            #set up some inital conditions and the values of the driving term f

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

            simul = SEM_Wave_1d.SEM_Wave(nodes, nsteps, Tend, N, fVals, omega)
            valueMatrix = SEM_Wave_1d.LaplaceMMS(simul, 1e-6, "polynomial") 

            #diffs[j] = maximum(abs.(valueMatrix[1] - valueMatrix[2]))
            index = Int(floor(length(valueMatrix)*0.6))
            #diffs[j] = maximum(abs.(valueMatrix[1] - valueMatrix[2]))
            diffs[j] = maximum(abs.(valueMatrix[1, index] - valueMatrix[2, index]))

        end

        plot!(numberOfNodesList, diffs, xscale=:log10, yscale=:log10, xlims=(numberOfNodesList[1], numberOfNodesList[end]), label = "N = " * string(N), legend=:bottomleft)            
        #plot!(numberOfNodesList[1:end-1], diffs, xscale=:log10, yscale=:log10, xlims=(numberOfNodesList[1], numberOfNodesList[end]), label = "N = " * string(N), legend=:bottomright)            
        println("done with N = " * string(N))        

    end

    savefig(plt, "mmsplot.png")

end



