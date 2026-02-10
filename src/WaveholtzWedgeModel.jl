


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
    xr = 600.0
    yl = 0.0
    yr = 1000.0



    N = 5   # degree of interpolation polynomials
    Kx = 90 # number of elements in x-direction
    Ky = Integer(ceil(Kx* (5/3))) # number of elements in y-direction

    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.
    #c_square(x, y) = 2100

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    
    omega = 40.0*pi
    a = 1/sqrt(2)
    b = sqrt(1-a^2)
    bc = [a, b]

    delta(x) = Float64(x==0.0) # the best implementation known to man
    f(x, y) = omega^2 * delta(abs(x - 300.0))' * delta(abs(y))
    fVals = f.(simul.x', simul.y)
    fVals = omega^2 * exp.(-0.1*(simul.y).^2) * exp.(-0.1*(simul.x .- 300.0).^2)'
    #fVals = omega^2 * exp.(-0.01*(simul.y .- 0.5*(yr+yl)).^2) * exp.(-0.01*(simul.x .- 0.5*(xr+xl)).^2)'
    

    g = zeros(length(simul.y), length(simul.x))

    u_0, u_1, nIters = SEM_Wave_2d.WaveholtzGMRES(simul, omega, -fVals, bc, g, 1e-6)
    #SEM_Wave_2d.Waveholtz(simul, omega, fVals, bc, g, 0.001)
    #SEM_Wave_2d.WaveholtzAnimation(simul, omega, fVals, bc, g, 30, true)
    #u_WHI = simul.uFiltered

    
    plt = heatmap(simul.x, simul.y, log10.(abs.(u_0)), aspect_ratio =:1.00)
    #plt = heatmap(simul.x, simul.y, log10.(abs.(u_WHI)))
    
    #plt = surface(simul.x, simul.y[end:-1:1], simul.c_square)
    #plt = heatmap(simul.x, simul.y, simul.c_square)
    #plt = heatmap(simul.x, simul.y, fVals)
    savefig(plt, "WaveholtzWedge")
    

    #lapU = SEM_Wave_2d.LaplaceTerm(simul, simul.uFiltered) / simul.timestep^2 + simul.M_b .* simul.uFiltered
    #diff = abs.(lapU + omega^2*simul.uFiltered - fVals)

    #tester = sin.(pi .* simul.x').*cos.(pi .* simul.y)
    #testerLap = SEM_Wave_2d.LaplaceTerm(simul, tester) / simul.timestep^2 + simul.M_b .* tester

    #relDiff = SEM_Wave_2d.LpNorm(simul, lapU + omega^2*simul.uFiltered, fVals, 2)/SEM_Wave_2d.LpNorm(simul, simul.uFiltered, zeros(length(simul.y), length(simul.x)), 2)
    #println(maximum(diff))
    #println(relDiff)

    #plt = heatmap(simul.x, simul.y, log10.(diff))
    #plt = heatmap(simul.x, simul.y, log10.(diff))
    #plt = heatmap(simul.x, simul.y, log10.(abs.(testerLap + pi^2 .* tester)))
    #plt = surface(simul.x, simul.y[end:-1:1], testerLap + pi^2 .* tester)
    #plt = surface(simul.x, simul.y[end:-1:1], diff)
    #savefig(plt, "WaveholtzDiff")

end

