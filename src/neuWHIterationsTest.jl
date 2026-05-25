


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

    N = 4
    
    tol = 1e-6
    omega = 9.225806451612902

    c_square(x, y) = 1
    c_square(x, y) = 1 - 0.9*exp.(-((y .- 0.1)/0.1).^2) * exp.(-((x .- 0.7)/0.1).^2)' - 0.9*exp.(-((y .- 0.3)/0.1).^2) * exp.(-((x .- 0.27)/0.1).^2)'

    c_min = 0.1 ### has to be set properly! ###

    L = maximum([xr-xl, yr-yl]) # length of the domain
    
    # rule of thumb for 4th order, rounded up a bit
    dofs_x = Int(ceil(1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
    numberOfElements = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!


    Kx = numberOfElements # number of elements in x-direction
    Ky = numberOfElements # number of elements in y-direction

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    #fVals = 10*exp.(-((simul.y .- 0.5)/0.01).^2) * exp.(-((simul.x .- 0.5)/0.01).^2)'
    fVals = cos.(2*pi*simul.y) * cos.(2*pi*simul.x)'

    g = zeros(length(simul.y'), length(simul.x))        
 
    #sol1, nIter = SEM_Wave_2d.Waveholtz(simul, omega, fVals, bc, g, tol)
    #u_0, u_1, history = SEM_Wave_2d.WaveholtzGMRESnew(simul, omega, fVals, bc, g, 1e-9)
    SEM_Wave_2d.WaveholtzAnimation(simul, omega, fVals, bc, g, 300)



    #SEM_Wave_2d.WaveholtzConvHistory(simul, omega, fVals, bc, g, 500)

    #println(dofs_x)
    #println(numberOfElements)

    #plt = surface(simul.y, simul.x[end:-1:1], u_0)
    #plt = heatmap(simul.x, simul.y, log10.(abs.(sol1)))
    #savefig(plt, "neuWHIterationsTest")

end







function main2()

    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    #a = 1/sqrt(2)
    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    N = 4
    
    tol = 1e-2
    omega = 9.225806451612902

    c_square(x, y) = 1
    c_min = 1.0 ### has to be set properly! ###

    L = maximum([xr-xl, yr-yl]) # length of the domain
    
    # rule of thumb for 4th order, rounded up a bit
    dofs_x = Int(ceil(1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
    numberOfElements = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!


    Kx = numberOfElements # number of elements in x-direction
    Ky = numberOfElements # number of elements in y-direction

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    fVals = 10*exp.(-((simul.y .- 0.5)/0.01).^2) * exp.(-((simul.x .- 0.5)/0.01).^2)'
    g = zeros(length(simul.y), length(simul.x))        
 


    ### from SEM_Wave_2d.Waveholtz:
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    spaceStep = minimum([delta_x delta_y])
    nsteps = Integer(ceil(1.5*Tend * (1/delta_x + 1/delta_y)))


    # do a lil' waveholtzin' and measure errors

    u_0, u_1, history = SEM_Wave_2d.WaveholtzGMRESnew(simul, omega, fVals, bc, g, 1e-9)

    k = 200
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))

    anim = Animation()

    for j = 1:k


        useMMS = false
        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)
        
        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

        println(maximum(log10.(abs.(simul.uFiltered - u_0))))
        #surface(simul.x, simul.y[end:-1:1], log10.(abs.(simul.uFiltered - u_0)), zlims=(-15, 1), legend=:false)
        surface(simul.x, simul.y[end:-1:1], abs.(simul.uFiltered - u_0), legend=:false)
        frame(anim)

        println(string(j) * " out of " * string(k))
    
    end

    gif(anim, "neuWHIterationsTest.gif", fps=10)

end

