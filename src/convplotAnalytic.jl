
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

    numberOfRuns = 10
    numberOfNodesList = zeros(numberOfRuns)
    #timestepsList = zeros(numerOfRuns)
    timestepsBase = 20

    for j = 1:numberOfRuns
        numberOfNodesList[j] = 2^j + 1
        #timestepsList[j] = 40*2^j + 1
    end

    numberOfNodesList = Int.(numberOfNodesList)
    #timestepsList = Int.(timestepsList)

    println(numberOfNodesList)

    ordersList = [3, 4, 5]#, 6, 7]
    
    omega = 6.0
    omega = 5.0*pi

    #GenerateConvPlot(ordersList, numberOfNodesList, timestepsList, numberOfRuns, omega)
    GenerateConvPlot(ordersList, numberOfNodesList, timestepsBase, numberOfRuns, omega)
    
end


function GenerateConvPlot(ordersList, numberOfNodesList, timestepsBase, numberOfRuns, omega)

    plt = plot()

    for n = 1:length(ordersList)
        #the number of quadrature points per element; frequency parameter
        N = ordersList[n]
        
        #number of timesteps + outer limit of time interval
        Tend = 2.0
        nsteps = timestepsBase
        #nsteps = timestepsList[n]

        #outer limits of the interval we are solving the equation on
        xl = 0
        xr = 1

        #where we store the difference between the jth run and the final one
        diffs = zeros(numberOfRuns)

        for j = 1:numberOfRuns
            #solve the wave equation with the different numbers of nodes
            
            #set up some inital conditions and the values of the driving term f

            #newSteps = Int(ceil((N/(N+1))*timestepsList[j]))
            #factor = 2^((N-2)/2)

            factor = 2
            newSteps = Int(ceil(nsteps*factor))
            nsteps = newSteps
            
            #nsteps = timestepsList[j]

            numberOfNodes = numberOfNodesList[j]
            
            nodes = zeros(numberOfNodes)
            
            range = LinRange(xl, xr, numberOfNodes)
            
            for i = 1:numberOfNodes
                nodes[i] = range[i]
            end

            quadpoints, ~ = gausslobatto(N)
            totalPointCount = (length(nodes)-1) * (N-1) + 1


            p = 3
            q = 3
            
            

            x = SEM_Wave_1d.ConstructX(nodes, quadpoints)
            fVals = RightHandSide(x, p, q, omega)
            #fVals = RightHandSide(x, 5, omega)
            


            uStart = Reference(x, 0, p, q, omega)
            #uStart = Reference2(x, 0, 5, omega)
            #uStartDer = StartDer(x, p, q, omega) this derivative is always zero
            uStartDer = zeros(totalPointCount)
            animate = false
            snapshotFrequency = 0   #as animate = false, this value is not used anywhere

            println(nsteps)
            println(length(nodes))
            output = SEM_Wave_1d.Simulate(nodes, nsteps, Tend, N, fVals, omega, uStart, uStartDer, animate, snapshotFrequency)

            diffs[j] = maximum(abs.(output - Reference(x, Tend, p, q, omega)))
            #diffs[j] = maximum(abs.(output - Reference2(x, Tend, 5, omega)))
            

  #=          
            if (N == ordersList[end] && j == numberOfRuns)
                plot!(x, output, xlims=(x[1], x[end]), label = "ouput")            
                plot!(x, Reference(x, Tend, p, q, omega), xlims=(x[1], x[end]), label = "reference")                            
                #plot!(x, Reference2(x, Tend, 5, omega), xlims=(x[1], x[end]), label = "reference")

            end
            
=#

        end

        println(diffs)

        plot!(numberOfNodesList, diffs, xscale=:log10, yscale=:log10, xlims=(numberOfNodesList[1], numberOfNodesList[end]), label = "N = " * string(N))            
        println("done with N = " * string(N))        

    end

    savefig(plt, "convplotAnalytic")
    #display(plt)

end


function Reference(x, t, p, q, omega)
    #outputs the reference solution of the type x^p*(1-x)^q * cos(\omega*t) evaluated at x, t
    out = (x.^p .* (-x .+ 1).^q) .* cos(omega*t)
end

function Reference2(x, t, n, omega)
    #outputs the reference solution of the type x^p*(1-x)^q * cos(\omega*t) evaluated at x, t
    out = 10*sin.(n*pi*x) .* cos(omega*t)
    #out = 1
end

function RightHandSide(x, p, q, omega)
    #the Refereence_xx - Reference_tt, evaluated at x, t
    out = (p*(p-1)*x.^(p-2).*(-x .+ 1).^q - 2*p*q*x.^(p-1).*(-x .+ 1).^(q-1) + q*(q-1)*x.^p .*(-x .+ 1).^(q-2)) + 
           omega^2*(x.^p .* (-x .+ 1).^q)
    #out = 1
           #SEM_Wave_1d already multiplies all this by cos(omega*t)
end

function RightHandSide2(x, n, omega)
    #the Refereence_xx - Reference_tt, evaluated at x, t
    out = (omega^2 - pi^2*n^2)*sin.(n*pi.*x)
    #out = 1
           #SEM_Wave_1d already multiplies all this by cos(omega*t)
end


function StartDer(x, p, q, omega)
    #gives the time derivative of Reference at t = 0
    out = -omega*(x.^p .* (-x .+ 1).^q) .* sin(omega*t)
    #out = 1
end
