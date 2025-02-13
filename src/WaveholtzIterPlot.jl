

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra

include("SEM_Wave_1d_new.jl")
using .SEM_Wave_1d


include("MMS.jl")
using .MMS




function main()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    xl = 0
    xr = 1

    #delta = 0.2
    #bc = [1/sqrt(2) + delta, sqrt(1/2 - sqrt(2)*delta - delta^2)]

    bc = [0.0, 1.0]
    #bc = [1/sqrt(2), 1/sqrt(2)]
    
    tol = 1e-9

    #set up range of omegas
    omegaSmall = 0.1
    omegaBig = 10.0
    numberOfOmegas = 20
    omegas = collect(LinRange(omegaSmall, omegaBig, numberOfOmegas))
    
    N = 3 # number of points we interpolate in for each element
    
    
    #run for these values of omegas, record the number of iterations
    nIterWaveholtz = zeros(numberOfOmegas)
    nIterGMRES = zeros(numberOfOmegas)
    
    for j = 1:numberOfOmegas

        omega = omegas[j]
        numberOfNodes = Int(ceil(10 + omega^(1 + 1/(N+1))))
        

        nodes = collect(LinRange(xl, xr, numberOfNodes))
        simul = SEM_Wave_1d.SEM_Wave(nodes, N)

        ###### set up the problem ######
        fVals = exp.(-50*(simul.x .- (xr-xl)/3).^2)

        simul.g = [0.0, 0.0]
        simul.bc = bc


        println("Waveholtz: " * string(omega))
        sol1, nIter1 = SEM_Wave_1d.Waveholtz(simul, omega, fVals, tol)

        println("Waveholtz + GMRES: " * string(omega))
        sol2, nIter2 = SEM_Wave_1d.WaveholtzGMRES(simul, omega, fVals, tol)

        nIterWaveholtz[j] = nIter1
        nIterGMRES[j] = nIter2
    
    end

    
    # plot the result

    plt = plot(omegas, nIterWaveholtz, yscale=:log10, label="Waveholtz", legend=:bottomleft)
    plot!(omegas, nIterGMRES, yscale=:log10, label="Waveholtz with GMRES", legend=:bottomleft, xaxis=:"ω", yaxis=:"number of iterations")
    
    savefig(plt, "WaveholtzIterPlot")

end

