
using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using TimerOutputs


include("SEM_Wave_1d.jl")





#=
Here we conduct an MMS-test for the Laplacian
=#
function main()

    #to = TimerOutput()


    ##########################
    #set up the time interval, the number of elements, the degree of the interpolation polynomials

    Tend = 10.0
    nsteps = 10000

    xl = 0
    xr = 1
    NumberOfNodes = 3200
    nodes = zeros(NumberOfNodes)
    
    range = LinRange(xl, xr, NumberOfNodes)
    
    for i = 1:NumberOfNodes
        nodes[i] = range[i]
    end

    N = 6 #number of points we interpolate in for each element
    

    #set up some inital conditions and the values of the driving term f
    totalPointCount = (length(nodes)-1) * (N-1) + 1
    fVals = zeros(totalPointCount)    
    omega = 2.0

    simul = SEM_Wave_1d.SEM_Wave(nodes, nsteps, Tend, N, fVals, omega)



    #@timeit to "Laplacian Computation" valueMatrix = SEM_Wave_1d.LaplaceMMS(simul, 1e-8, "polynomial") 

    valueMatrix = SEM_Wave_1d.LaplaceMMS(simul, 1e-8, "polynomial") 

    plt = plot(simul.x[2:end-1], [valueMatrix[1] valueMatrix[2]], xlims=(simul.x[2], simul.x[end-1]), ylims=(-maximum(maximum(valueMatrix)), maximum(maximum(valueMatrix))), label = ["using LaplaceCalculation" "computed symbolically"])
    #plt = plot(x[2:end-1], abs.(valueMatrix[1] - valueMatrix[2]), xlims=(x[2],x[end-1]))#, ylims=(-10, 15), label = ["using LaplaceCalculation" "computed symbolically"])
    #plt = plot(x[2:end-1], [valueMatrix[1] valueMatrix[2]], xlims=(x[2],x[end-1]))
    plot!(simul.x[2:end-1], zeros(length(simul.x[2:end-1])), seriestype=:scatter, ms=0.5, label = "quadrature points")

    savefig(plt, "testplot.png")



    println(maximum(abs.(valueMatrix[1] - valueMatrix[2])))

    #show(to)

end
