


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
    


    N = 8   # degree of interpolation polynomials
    Kx = 30 # number of elements in x-direction
    Ky = 30 # number of elements in y-direction

    heaviside(x) = 0.5 * (sign(x) + 1)

    #c_square(x, y) = 1 - 0.999*((heaviside(x-0.2) - heaviside(x-0.3)) .* (heaviside(y+0.9) - heaviside(y-0.9)) + (heaviside(x+0.3) - heaviside(x+0.2)) .* (heaviside(y+0.9) - heaviside(y-0.9)))
    #c_square(x, y) = 1 - 0.9999*((heaviside(x-0.025) - heaviside(x-1)) .* (heaviside(y+0.5) - heaviside(y)) + (heaviside(x+1) - heaviside(x+0.025)) .* (heaviside(y+0.5) - heaviside(y)))
    #c_square(x, y) = (1 - 0.999*((heaviside(x-0.305) - heaviside(x-1)) .* (heaviside(y+0.5) - heaviside(y)) + (heaviside(x+0.295) - heaviside(x-0.295)) .* (heaviside(y+0.5) - heaviside(y)) + (heaviside(x+1) - heaviside(x+0.305)) .* (heaviside(y+0.5) - heaviside(y)))).^2
    c_square(x, y) = 1 - 0.9*exp.(-((y .- 0.1)/0.1).^2) * exp.(-((x .+ 0.7)/0.1).^2)' - 0.9*exp.(-((y .+ 0.3)/0.1).^2) * exp.(-((x .- 0.27)/0.1).^2)'
    #c_square(x, y) = 2 - 1.9*exp.(-(y/0.05).^2) * exp.(-(x/0.1).^2)' - 1.9*exp.(-((y .- 1.0)/0.05).^2) * exp.(-((x .- 1.0)/0.1).^2)'
    #c_square(x, y) = 1.0


    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    fVals = 10*exp.(-((simul.y .- (yr+yl)/2)/0.1).^2) * exp.(-((simul.x .- (xr+xl)/2)/0.1).^2)'
    #fVals = 10*exp.(-((simul.y .- 0.3)/0.01).^2) * exp.(-((simul.x .+ 0.6)/0.01).^2)' + 10*exp.(-((simul.y .+ 0.3)/0.01).^2) * exp.(-((simul.x .+ 0.6)/0.01).^2)'
    #fVals = 10*exp.(-((simul.y .- 0.5)/0.01).^2) * exp.(-((simul.x .- 0.5)/0.01).^2)'
    #fVals = zeros(length(simul.y), length(simul.x))
    
    omega = 4.0*pi

    
    #fVals = omega^2*exp.(-((simul.y .+ 0.8)/0.01).^2) * exp.(-((simul.x .- 0.3)/0.01).^2)' + omega^2*exp.(-((simul.y .+ 0.8)/0.01).^2) * exp.(-((simul.x .+ 0.3)/0.01).^2)'

    #epsilon = 1e-2
    #omega = 2*pi + epsilon
    
    a = 1/sqrt(2)
    #a = 0.0
    b = sqrt(1-a^2)
    bc = [a, b]

    g = zeros(length(simul.y), length(simul.x))

    
    #SEM_Wave_2d.Waveholtz(simul, omega, fVals, bc, g, 4)
    #SEM_Wave_2d.Waveholtz(simul, omega, fVals, bc, g, 0.001)
    #u_WHI = simul.uFiltered
    
    #u_0, u_1, nIters = SEM_Wave_2d.WaveholtzGMRES(simul, omega, fVals, bc, g, 1e-6)
    
    nIter = 40
    u_0, u_1, history = SEM_Wave_2d.WaveholtzConvGMRES(simul, omega, fVals, bc, g, nIter)
    println(history)

    
    #SEM_Wave_2d.WaveholtzAnimation(simul, omega, fVals, bc, g, 50)
    #SEM_Wave_2d.WaveholtzConvHistory(simul, omega, fVals, bc, g, 500)
    
    #res1 = SEM_Wave_2d.ErrorEstimate(simul, u_WHI, fVals, omega)
    #res2 = SEM_Wave_2d.ErrorEstimate(simul, u_0, fVals, omega)

    #println(res1)
    #println(res2)
    #println(norm(u_WHI - u_0))

    #plt = surface(simul.x, simul.y[end:-1:1], log10.(abs.(simul.uFiltered)))
    #plt = surface(simul.x, simul.y[end:-1:1], simul.uFiltered)
    #plt = heatmap(simul.x, simul.y, log10.(abs.(simul.uFiltered)))
    plt = heatmap(simul.x, simul.y, log10.(abs.(u_0)))
    #plt = heatmap(simul.x, simul.y, log10.(abs.(u_WHI)))
    
    #plt = surface(simul.x, simul.y[end:-1:1], simul.c_square)
    #plt = surface(simul.x, simul.y[end:-1:1], fVals)
    savefig(plt, "WaveholtzTest")

end

