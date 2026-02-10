

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_2d.jl")
using .SEM_Wave_2d


include("MMS.jl")
using .MMS




function main()

    
    xl = -2.123/5
    xr = 4.0234/5
    yl = -3.234567/5
    yr = 1.5678/5
    

    #=
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0
    =#


    numberOfRuns = 5
    numbersOfElements = zeros(numberOfRuns, 1)
    diffs = zeros(numberOfRuns, 1)

    Ns = [3 4 5]
    plt = plot()

    for i = length(Ns):-1:1
        
        N = Ns[i]

        for j = 1:numberOfRuns
            
            #Kx = Integer(6 + 2*j)
            Kx = Integer(round(2^((j+1)/2 + 2))) # number of elements in x-direction

            Ky = Kx # number of elements in y-direction
            numbersOfElements[j] = Kx
            println(Kx)

            #bc = [0.0, 1.0]
            #bc = [1/sqrt(2), 1/sqrt(2)]
            a = 0.1234
            #a = 0.0
            b = sqrt(1-a^2)
            bc = [a, b]

            c_square(x, y) = 1 

            simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

            g = zeros(length(simul.y), length(simul.x))
            forcing = zeros(length(simul.y), length(simul.x))
            Tend = 1.0
            nsteps = 200000
            omega = 1.234
            


            simul.useMMS = true
            simul.MMS_j.type = [3 3 2]

            simul.MMS_j.coeff[1, 1] = 1
            simul.MMS_j.coeff[1, 2] = 3*pi * 2 / (xr - xl)
            #simul.MMS_j.coeff[1, 2] = 0.0
            #simul.MMS_j.coeff[1, 3] = -3*pi*(xr + xl) / (xr - xl)
            simul.MMS_j.coeff[1, 3] = 5.0
            #simul.MMS_j.coeff[1, 3] = 6.0
            
            simul.MMS_j.coeff[2, 1] = 1
            #simul.MMS_j.coeff[2, 2] = 5*pi * 2 / (yr - yl)
            simul.MMS_j.coeff[2, 2] = 0.5
            #simul.MMS_j.coeff[2, 3] = 0.0
            #simul.MMS_j.coeff[2, 3] = -5*pi*(yr + yl) / (yr - yl)
            simul.MMS_j.coeff[2, 3] = 2
            #simul.MMS_j.coeff[2, 3] = 6.0

            simul.MMS_j.coeff[3, 1] = 1
            #simul.MMS_j.coeff[3, 2] = 2*pi
            simul.MMS_j.coeff[3, 2] = 1
            #simul.MMS_j.coeff[3, 3] = 0
            simul.MMS_j.coeff[3, 3] = 2

            simul.MMS_c.type = [3 3]

            simul.MMS_c.coeff[1, 1] = 1
            simul.MMS_c.coeff[1, 2] = 0
            simul.MMS_c.coeff[1, 3] = 1
            #simul.MMS_c.coeff[1, 1] = 1
            simul.MMS_c.coeff[1, 2] = (xr+xl)/2
            #simul.MMS_c.coeff[1, 3] = 0.05


            simul.MMS_c.coeff[2, 1] = 1
            #simul.MMS_c.coeff[2, 2] = yl-1
            simul.MMS_c.coeff[2, 2] = 0.0
            simul.MMS_c.coeff[2, 3] = 2
            #simul.MMS_c.coeff[2, 1] = 1
            #simul.MMS_c.coeff[2, 2] = (yr+yl)/2
            #simul.MMS_c.coeff[2, 3] = 0.05




            #=
            simul.useMMS = true
            simul.MMS_j.type = [3 3 3]

            simul.MMS_j.coeff[1, 1] = 1
            simul.MMS_j.coeff[1, 2] = 0.0
            simul.MMS_j.coeff[1, 3] = 5.0
            
            simul.MMS_j.coeff[2, 1] = 1
            simul.MMS_j.coeff[2, 2] = 0.0
            simul.MMS_j.coeff[2, 3] = 0.0

            simul.MMS_j.coeff[3, 1] = 1
            simul.MMS_j.coeff[3, 2] = 1
            simul.MMS_j.coeff[3, 3] = 0

            simul.MMS_c.type = [3 3]

            simul.MMS_c.coeff[1, 1] = 1
            simul.MMS_c.coeff[1, 2] = 0
            simul.MMS_c.coeff[1, 3] = 1

            simul.MMS_c.coeff[2, 1] = 1
            simul.MMS_c.coeff[2, 2] = 0.0
            simul.MMS_c.coeff[2, 3] = 2

            =#


            simul.MMS_j.type = [2 2 2]
            simul.MMS_j.coeff[1, 1] = 1.98765
            simul.MMS_j.coeff[1, 2] = 1 + 0.57721566490153286060651209008240243104215933593992
            simul.MMS_j.coeff[1, 3] = log(20)

            
            
            simul.MMS_j.coeff[2, 1] = 0.12345
            simul.MMS_j.coeff[2, 2] = sqrt(17)
            simul.MMS_j.coeff[2, 3] = 10*exp(-pi)

            

            simul.MMS_j.coeff[3, 1] = (1+sqrt(5))/2
            simul.MMS_j.coeff[3, 2] = 5*pi/4
            simul.MMS_j.coeff[3, 3] = 4


            simul.MMS_c.type = [3 3]

            simul.MMS_c.coeff[1, 1] = 1
            simul.MMS_c.coeff[1, 2] = (xr+xl)/2
            simul.MMS_c.coeff[1, 3] = 2
            


            simul.MMS_c.coeff[2, 1] = 1
            simul.MMS_c.coeff[2, 2] = (yr+yl)/2
            simul.MMS_c.coeff[2, 3] = 2





            # for use in tests of the initialisation
            uStart = zeros(length(simul.y), length(simul.x))
            uStartDer = zeros(length(simul.y), length(simul.x))

            for s = 1:length(simul.y)
                for t = 1:length(simul.x)
                    uStart[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 0, simul.MMS_j)
                    uStartDer[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 1, simul.MMS_j)
                end
            end

            #mmsRef, appx = SEM_Wave_2d.ForcingMMS(simul, 10) # stepnumber 10 chosen arbitrarily.
            #mmsRef, appx = SEM_Wave_2d.LaplaceMMS(simul)
            mmsRef, appx = SEM_Wave_2d.InitialiseMMS(simul, uStart, uStartDer, Tend, nsteps, forcing, omega, bc, g)
            #mmsRef, appx = SEM_Wave_2d.BoundaryMMS(simul, 10, bc)

            #println(maximum(mmsRef))
            #println(maximum(appx))
            diffs[j] = maximum(abs.(mmsRef[2:end-1, 2:end-1]-appx[2:end-1, 2:end-1]))
            #diffs[j] = SEM_Wave_2d.LpNorm(simul, mmsRef, appx, 2) / SEM_Wave_2d.LpNorm(simul, mmsRef, zeros(length(simul.y), length(simul.x)), 2) 
            println(maximum(abs.(mmsRef-appx)))

            #plt2 = surface(simul.x[2:end-1], simul.y[end-1:-1:2], log10.(abs.((mmsRef-appx)[2:end-1, 2:end-1]./mmsRef[2:end-1, 2:end-1])))
            #savefig(plt2, "Lap_tester")

            

        end

        scatter!(numbersOfElements, diffs, xscale=:log10, yscale=:log10, label = "N = " * string(N), ms = 3)

        logLine = zeros(numberOfRuns, 1)
        logLine[1] = diffs[1]

        referenceOrder = N

        for j = 2:numberOfRuns
            logLine[j] = logLine[j-1] * (numbersOfElements[j-1] / numbersOfElements[j])^(referenceOrder)
        end

        plot!(numbersOfElements, logLine, xscale=:log10, yscale=:log10, label = "reference, order = " * string(referenceOrder), linestyle=:dash)


        println(diffs)

        #plot!(numbersOfElements, logLine, xscale=:log10, yscale=:log10, label = "reference, order = " * string(N+1), linestyle=:dash)

    end

    
    
    savefig(plt, "mmsLap_vs_appxLap")
    


#=
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    
    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    
    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, true, 1, 10)
=#

    

end

