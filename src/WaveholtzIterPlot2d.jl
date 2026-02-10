


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
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    #a = 1/sqrt(2)
    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    
    tol = 1e-3

    #set up range of omegas
    omegaSmall = 2.0
    omegaBig = 10.0
    numberOfOmegas = 40
    omegas = collect(LinRange(omegaSmall, omegaBig, numberOfOmegas))

    N = 5 # number of points we interpolate in for each element
    
    #run for these values of omegas, record the number of iterations
    nIterWaveholtz = zeros(numberOfOmegas)


    for j = 1:numberOfOmegas

        omega = omegas[j]
        #numberOfNodes = Int(ceil(10 + omega^(1 + 1/(N+1))))
        numberOfNodes = 10
        

        Kx = numberOfNodes # number of elements in x-direction
        Ky = numberOfNodes # number of elements in y-direction

        #c_square(x, y) = 1
        c_square(x, y) = 1 - 0.9*exp.(-((y .- 0.15)/0.05).^2) * exp.(-((x .- 0.7)/0.1).^2)' - 0.9*exp.(-((y .- 0.8)/0.05).^2) * exp.(-((x .- 0.27)/0.1).^2)'


        simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

        ###### set up the problem ######
        #fVals = exp.(-((simul.y .- (yr+yl))/0.1).^2) * exp.(-((simul.x .- (xr+xl))/0.1).^2)'
        #fVals = 10*exp.(-((simul.y .- 0.3)/0.01).^2) * exp.(-((simul.x .+ 0.6)/0.01).^2)' + 10*exp.(-((simul.y .+ 0.3)/0.01).^2) * exp.(-((simul.x .+ 0.6)/0.01).^2)'
        fVals = 10*exp.(-((simul.y .- 0.5)/0.01).^2) * exp.(-((simul.x .- 0.5)/0.01).^2)'
        g = zeros(length(simul.y), length(simul.x))        

        println("omega: " * string(omega))
        sol1, nIter = SEM_Wave_2d.Waveholtz(simul, omega, fVals, bc, g, tol)

        nIterWaveholtz[j] = nIter
    
    end

    
    # plot the result

    plt = plot(omegas, nIterWaveholtz, yscale=:log10, label="Waveholtz", legend=:bottomleft)
    
    savefig(plt, "WaveholtzIterPlot2d")

    println("values of omega: ")
    println(omegas)
    println("iterations for convergence: ")
    println(nIterWaveholtz)

end

