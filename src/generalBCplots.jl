

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

    xl = 0
    xr = 1
    numberOfLevels = 5
    numbersOfNodes = zeros(numberOfLevels)

    g = [0.0, 0.0]
    omega = 8.13

    delta = 0.1
    bc = [1/sqrt(2) + delta, sqrt(1/2 - sqrt(2)*delta - delta^2)]    

    #bc = [0.0, 1.0]
    #bc = [1/sqrt(2), 1/sqrt(2)]
    
    for i = 1:numberOfLevels

        numbersOfNodes[i] = Integer(5*2^i + 1)  #the number of elements is 10*2^i
        
    end

    numbersOfNodes = Integer.(numbersOfNodes)


    Ns = [3 4 5] #number of points we interpolate in for each element


    plt = plot()
    reference = zeros(numbersOfNodes[end])
    

    
    for N in Ns

        diffs = zeros(numberOfLevels - 1)
        #diffs2 = zeros(numberOfLevels - 1)
        #x_ref = zeros((numbersOfNodes[end]-1)*N + 1)

        for j = length(numbersOfNodes):-1:1
            
            nodes = collect(LinRange(xl, xr, numbersOfNodes[j]))

            simul = SEM_Wave_1d.SEM_Wave(nodes, N)

            simul.Tend = 2.47
            #simul.Tend = 0.22
            #simul.nsteps = 10*numbersOfNodes[j]*N
                        
            delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
            simul.nsteps = Integer(ceil(50*simul.Tend/delta_x))
            simul.timestep = simul.Tend/simul.nsteps
            
            #println("CFL number: " * string(simul.timestep/delta_x))
            

            simul.omega = omega
            simul.g = g
            simul.bc = bc
            
            #simul.fVals = zeros(length(simul.x))
            simul.fVals = 150*exp.(-((simul.x.-0.5)/0.1).^2)

            uStart = 5*exp.(-((simul.x.-0.5)/0.1).^2)
            uStartDer = zeros(length(simul.x))

            animate = false
            
            
            SEM_Wave_1d.Simulate(simul, uStart, uStartDer, animate)


            ### Compare to finest run ###
            if (j == length(numbersOfNodes))
                reference = simul.uNow
                #x_ref = simul.x
            else    

                #println(maximum(simul.x[1:N-1:end] -x_ref[1:(N-1)*2^(length(numbersOfNodes) - j):end]))

                diffs[j] = maximum(abs.(reference[1:(N-1)*2^(length(numbersOfNodes) - j):end] - simul.uNow[1:N-1:end]))
                #diffs2[j] = maximum(abs.(simul.uNow[1:N-1:end]))

            end

        end

        println("order: " * string(N))
        
        for i = 1:length(diffs)-1
            #println(log(diffs[i]/diffs[i+1])/log(2))
            println(diffs[i]/diffs[i+1])
        end

        #println(diffs)
        plt = plot!(numbersOfNodes[1:end-1], diffs, xscale=:log10, yscale=:log10, label = "N = " * string(N))
        #plt = plot!(numbersOfNodes[1:end-1], diffs2, xscale=:log10, yscale=:log10, label = "test: N = " * string(N))

    end


    savefig(plt, "convPlotsGeneralBC")


end

