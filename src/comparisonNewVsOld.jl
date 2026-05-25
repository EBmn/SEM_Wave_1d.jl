

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using Base.Threads

include("SEM_Wave_2d.jl")
using .SEM_Wave_2d


include("MMS.jl")
using .MMS



function mmsConvPlotOld()
    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    
    
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0
    

    N = 3 # degree of interpolation polynomials
    numberOfRuns = 5
    errorData = zeros(numberOfRuns, 1)
    KxVec = zeros(numberOfRuns, 1)
    
    alpha = 0.0
    beta = 1.0

    for j = numberOfRuns:-1:1

        c_square(x, y) = 1
        omega = 7.1
        Tend = 0.0234
        nsteps = Int(5*2^(j-1))
        #nsteps = 1
        println(nsteps)

        Kx = Int(2^(j)) # number of elements in x-direction
        Ky = Kx # number of elements in y-direction

        KxVec[j] = Kx

        println("$Kx by $Ky elements running")

    
        simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

        simul.useMMS = true # we mean to perform an mms test


        #=
        # polynomial test
        simul.MMS_j.type = [3 3 3]

        simul.MMS_j.coeff[1, 1] = 1
        #simul.MMS_j.coeff[1, 2] = 0.567
        simul.MMS_j.coeff[1, 2] = 0.0
        simul.MMS_j.coeff[1, 3] = 1.0

        simul.MMS_j.coeff[2, 1] = 1
        simul.MMS_j.coeff[2, 2] = 0.0
        simul.MMS_j.coeff[2, 3] = 0.0

        simul.MMS_j.coeff[3, 1] = 1
        simul.MMS_j.coeff[3, 2] = 0
        simul.MMS_j.coeff[3, 3] = 4
        =#


        
        # trigonometric test
        simul.MMS_j.type = [2 2 2]
        simul.MMS_j.coeff[1, 1] = 1.0
        simul.MMS_j.coeff[1, 2] = 2*pi
        simul.MMS_j.coeff[1, 3] = 0.0
        
        simul.MMS_j.coeff[2, 1] = 1
        simul.MMS_j.coeff[2, 2] = 2*pi
        simul.MMS_j.coeff[2, 3] = 0.0

        simul.MMS_j.coeff[3, 1] = 1.0
        simul.MMS_j.coeff[3, 2] = 2.0
        simul.MMS_j.coeff[3, 3] = 0.0
        


        ###### set up the problem ######
        fVals = 1e4*exp.(-((simul.y .- 0.2*(yr+yl))/0.03).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.03).^2)' # should not matter what we put here, as the MMS functionality should put the proper values for f regardless


        g = zeros(length(simul.y), length(simul.x))


        # construct initial data
        uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
        #uStart = exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'
        #uStart = 2 .- (abs.(simul.y) .+ abs.(simul.x'))

        uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
        #uStartDer = 1e3*exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'

        reference = zeros(length(simul.y), length(simul.x))

        if simul.useMMS
            for s = 1:length(simul.y)
                for t = 1:length(simul.x)
                    uStart[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 0, simul.MMS_j)
                    uStartDer[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 1, simul.MMS_j)
                    reference[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], Tend, 0, 0, 0, simul.MMS_j)
                end
            end
        end

        bc = [alpha, beta]
        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, false)

        #println("older: ")
        #println(simul.uNow)

        errorData[j] = maximum(abs.(simul.uNow - reference))

        #plt = surface(simul.x, simul.y[end:-1:1], abs.(simul.uNow - reference))    
        plt = surface(simul.x, simul.y[end:-1:1], simul.uNow - reference)    
        #plt = surface(simul.x, simul.y[end:-1:1], reference)
        savefig(plt, "testsurface" * string(j))

        println(maximum(reference))
        println(maximum(simul.uNow))


    end

    println(errorData)
    plt = plot(KxVec, errorData, xscale=:log10, yscale=:log10)
    savefig(plt, "mmsConvPlot")

end
