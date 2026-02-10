

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_2d.jl")
using .SEM_Wave_2d


include("MMS.jl")
using .MMS




function main()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0
    

    #delta = 0.25
    #bc = [1/sqrt(2) + delta, sqrt(1/2 - sqrt(2)*delta - delta^2)]

    bc = [0.0, 1.0]
    bc = [1/sqrt(2), 1/sqrt(2)]


    numberOfRuns = 4
    numbersOfElements = zeros(numberOfRuns)

    g = [0.0, 0.0]
    omega = 8.13

    
    for i = 1:numberOfRuns

        numbersOfElements[i] = Integer(2^(i+4))
        #numbersOfElements[i] = Integer(round(2^((i+1)/2+1)))
        
    end


    Ns = [3, 4]

    plt = plot()

    for i = length(Ns):-1:1

        N = Ns[i]  # degree of interpolation polynomials
        
        #c_square(x, y) = 1

        
        c_square(x, y) = 1 + (x-0.2).^2


        reference = zeros(Integer(numbersOfElements[end] * N + 1), (Integer(numbersOfElements[end] * N + 1)))

        diffs = zeros(numberOfRuns-1)

        for j = numberOfRuns:-1:1

            println(numbersOfElements[j])


            Kx = Integer(numbersOfElements[j]) # number of elements in x-direction
            Ky = Integer(numbersOfElements[j]) # number of elements in y-direction

            simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)


            #fVals = exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.1).^2)'
            fVals = zeros(length(simul.y), length(simul.x))
            g = zeros(length(simul.y), length(simul.x))

            omega = 7.1

            Tend = 0.25
            nsteps = 500

            uStart = exp.(-((simul.y .- 0.5*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.5*(xr+xl))/0.1).^2)'

            uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))

            SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, false)


            ### Compare to finest run ###

            if (j == length(numbersOfElements))
                reference = simul.uNow
            else    

                coarserReference = reference[1:N*2^(length(numbersOfElements) - j):end, 1:N*2^(length(numbersOfElements) - j):end]

                diffs[j] = maximum(abs.(simul.uNow[1:N:end, 1:N:end] - coarserReference))
                #diffs[j] = maximum(abs.(simul.uNow - reference))
                #diffs[j] = LpDist(simul.xNodes, simul.yNodes, simul.uNow[1:N:end, 1:N:end], coarserReference, 2)
                #diffs[j] = maximum(abs.(reference[1:N*2^(length(numbersOfElements) - j):end, 1:N*2^(length(numbersOfElements) - j):end] - simul.uNow[1:N:end, 1:N:end]))
                
            end

        end

        println(diffs)

        logLine = zeros(numberOfRuns-1, 1)
        
        #=
        logLine[end] = diffs[end]

        for j = numberOfRuns-2:-1:1
            logLine[j] = logLine[j+1] * (numbersOfElements[j+1] / numbersOfElements[j])^(N+1)
        end
        =#


        logLine[1] = diffs[1]

        for j = 2:numberOfRuns-1
            logLine[j] = logLine[j-1] * (numbersOfElements[j-1] / numbersOfElements[j])^(N+3)
        end


        println(numbersOfElements[1:end-1])
        scatter!(numbersOfElements[1:end-1], diffs, xscale=:log10, yscale=:log10, label = "N = " * string(N), ms = 3)
        plot!(numbersOfElements[1:end-1], logLine, xscale=:log10, yscale=:log10, label = "reference, order = " * string(N+1), linestyle=:dash)
    

    end
    
    savefig(plt, "convPlot2d")

end

function LpDist(xVals::Vector{Float64}, yVals::Vector{Float64}, u::Matrix{Float64}, v::Matrix{Float64}, p::Int64)

    # approximates the Lp-distance between u and v using the trapezoidal rule
    
    Kx = length(xVals) - 1                                       # number of points in the x-direction
    Ky = length(yVals) - 1                                       # number of points in the y-direction        
    
    delta_x = (xVals[end]-xVals[1])/Kx
    delta_y = (yVals[end]-yVals[1])/Ky
   
    weightsMat = (delta_x*delta_y) * [0.5; ones(Ky-1, 1); 0.5] * [0.5; ones(Kx-1, 1); 0.5]' 
    integrand = abs.(u-v).^p

    norm = sum(weightsMat .* integrand)^(1/p)

    return norm


end
