

using LinearAlgebra
using Plots
using FastGaussQuadrature

include("SEM_Wave_1d_new.jl")
using .SEM_Wave_1d


include("MMS.jl")
using .MMS



function main()

    ##########################
    #set up the time interval, the number of elements, the degree of the interpolation polynomials

    xl = 0
    xr = 1
    numberOfLevels = 8
    numbersOfNodes = zeros(numberOfLevels)

    #g = [0.0, 0.0]
    omega = 0.0
    

    delta = 0.2
    bc = [1/sqrt(2) + delta, sqrt(1/2 - sqrt(2)*delta - delta^2)]

    
    #bc = [0.0, 1.0]
    bc = [1/sqrt(2), 1/sqrt(2)]

    println(bc)
    
    for i = 1:numberOfLevels
        
        numbersOfNodes[i] = Integer(1*2^i + 1)
        
    end

    numbersOfNodes = Integer.(numbersOfNodes)


    Ns = [3, 4, 5] #number of points we interpolate in for each element


    plt = plot()
    reference = zeros(numbersOfNodes[end])



    for N in Ns

        diffs = zeros(numberOfLevels)
        diffs2 = zeros(numberOfLevels)
        lapDiff = zeros(numberOfLevels)
        fDiff = zeros(numberOfLevels)
        prevDiff = zeros(numberOfLevels)

        #used to draw a reference slope
        line = zeros(numberOfLevels)

        for j = 1:length(numbersOfNodes)

            nodes = collect(LinRange(xl, xr, numbersOfNodes[j]))

            simul = SEM_Wave_1d.SEM_Wave(nodes, N)

            simul.Tend = 2.47
            #simul.Tend = 0.22
            #simul.nsteps = 10*numbersOfNodes[j]*N

            delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
            simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
            simul.timestep = simul.Tend/simul.nsteps

            println(simul.nsteps)
            
            
            #println("CFL number: " * string(simul.timestep/delta_x))
            

            
            simul.bc = bc
            
            #set up MMS with a function that satisfies a Neumann condition
            simul.useMMS = true
            simul.MMS_j.type = 2
            

            simul.MMS_j.coeff[1, 2] = 13*pi/(simul.x[end] - simul.x[1])
            #simul.MMS_j.coeff[1, 2] = 12
            simul.MMS_j.coeff[1, 2] = 0
            simul.MMS_j.coeff[1, 3] = 0
            simul.MMS_j.coeff[2, 2] = omega
            simul.MMS_j.coeff[2, 3] = 0


            #simul.MMS_j.type = 3
            #simul.MMS_j.coeff[2, 1] = 4


            
            # the final boss of coefficients!
            #seems to work for the Neumann conditions
            simul.MMS_j.coeff[1, 2] = 13*pi/(simul.x[end] - simul.x[1])
            simul.MMS_j.coeff[1, 3] = 2.12
            simul.MMS_j.coeff[2, 2] = 1.17
            simul.MMS_j.coeff[2, 3] = -1.12
            


            #println(simul.MMS_j.coeff)


            uStart = zeros(length(simul.x))
            uStartDer = zeros(length(simul.x))
            reference = zeros(length(simul.x))
            
            for j = 1:length(simul.x)
                
                uStart[j] = MMS.MMSfun(simul.x[j], 0.0, 0, 0, simul.MMS_j)
                uStartDer[j] = MMS.MMSfun(simul.x[j], 0.0, 0, 1, simul.MMS_j)
                reference[j] = MMS.MMSfun(simul.x[j], simul.Tend, 0, 0, simul.MMS_j)

            end


            animate = false

            #compute the reference solution in the relevant points            
            SEM_Wave_1d.Simulate(simul, uStart, uStartDer, animate)


            ### Compare to reference in simul.MMS_j run ###            
            diffs[j] = maximum(abs.(reference - simul.uNow))
            #diffs2[j] = LpErr(simul, reference, 2)


            
            if j == length(numbersOfNodes)
                #plt = plot!(simul.x, reference, label = "ref N = " * string(N))
                #println(maximum(reference))
                #plt = plot!(simul.x, simul.uNow, label = "N = " * string(N))
            end

            if N == Ns[end]
                #plt = plot!(simul.x, log10.(abs.((simul.uNow - reference))), label = "N = " * string(N))
            end
            
        
            # Code for testing the Laplacian term using the MMS functionality: that bit seems to work
            #=
            lap, refLap = SEM_Wave_1d.LaplaceMMS(simul)
            #println(lap)
            lapDiff[j] = maximum(abs.(lap-refLap))
            order = N-1
            line[j] = lapDiff[1]*Float64(numbersOfNodes[1]^order)*Float64(numbersOfNodes[j])^(-order)
            

            if j == length(numbersOfNodes)
                #plt = plot!(simul.x[2:end-1], refLap, label = "ref N = " * string(N))
                #println(maximum(reference))
                #plt = plot!(simul.x[2:end-1], lap, label = "N = " * string(N))
            end
            =#

            # Code for testing the initialisation using the MMS functionality: really does seem to work...
            #=
            uPrev, uPrevRef = SEM_Wave_1d.InitialiseMMS(simul, uStart, uStartDer)
            
            prevDiff[j] = maximum(abs.(uPrev-uPrevRef))
            order = N-1
            line[j] = prevDiff[1]*Float64(numbersOfNodes[1]^order)*Float64(numbersOfNodes[j])^(-order)
            

            if j == length(numbersOfNodes)
                #plt = plot!(simul.x, uPrev, label = "ref N = " * string(N))
                #println(maximum(reference))
                #plt = plot!(simul.x, uPrevRef, label = "N = " * string(N))
            end
            =#


        end

        for j = 1:length(numbersOfNodes)
            
            if (mod(N, 2)==1)
                order = N+1
            else
                order = N+2
            end
            
            line[j] = diffs[5]*Float64(numbersOfNodes[5]^order)*Float64(numbersOfNodes[j])^(-order)

        end


        println("order: " * string(N))
        
        #=
        for i = 1:length(diffs)-1
            #println(log(diffs[i]/diffs[i+1])/log(2))
            println(diffs[i]/diffs[i+1])
        end
        =#
        #println(line)

        #println(diffs)
        plt = plot!(numbersOfNodes, diffs, xscale=:log10, yscale=:log10, label = "N = " * string(N))
        #plt = plot!(numbersOfNodes, diffs2, xscale=:log10, yscale=:log10, label = "N = " * string(N) * " L2", linestyle=:dot)
        plt = plot!(numbersOfNodes, line .+ 1e-15, xscale=:log10, yscale=:log10, label = false, linestyle=:dash, linecolor=:gray)

        



        #plt = plot!(numbersOfNodes, lapDiff, xscale=:log10, yscale=:log10, label = "laplace error: N = " * string(N))
        #plt = plot!(numbersOfNodes, fDiff .+ 1e-15, xscale=:log10, yscale=:log10, label = "forcing error: N = " * string(N))
        #plt = plot!(numbersOfNodes, prevDiff.+ 1e-15, xscale=:log10, yscale=:log10, label = "error in uPrev: N = " * string(N))
        #plt = plot!(numbersOfNodes, line, xscale=:log10, yscale=:log10, label = false, linestyle=:dash, linecolor=:gray)


    end

    savefig(plt, "mmsConvPlotsGeneralBC")

end
