


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

    a = 1/sqrt(2)
    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    N = 4
    
    tol = 1e-3
    omega = 7.0

    #c_square(x, y) = 1
    c_square(x, y) = 1 - 0.7*exp.(-((y .- 0.15)/0.05).^2) * exp.(-((x .- 0.7)/0.1).^2)' - 0.7*exp.(-((y .- 0.8)/0.05).^2) * exp.(-((x .- 0.27)/0.1).^2)'
    #c_square(x, y) = 1 - 0.9*exp.(-((y .- 0.15)/0.05).^2)

    c_min = 1.0 ### has to be set properly! ###

    L = maximum([xr-xl, yr-yl]) # length of the domain

    # rule of thumb for 4th order, rounded up a bit
    dofs_x = Int(ceil(2*1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
    numberOfElements = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!


    Kx = numberOfElements # number of elements in x-direction
    Ky = numberOfElements # number of elements in y-direction

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    
    fVals = 10*exp.(-((simul.y .- (yr+yl)/2)/0.1).^2) * exp.(-((simul.x .- (xr+xl)/2)/0.1).^2)'

    g = zeros(length(simul.y), length(simul.x))        

    ### from SEM_Wave_2d.Waveholtz:
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    spaceStep = minimum([delta_x delta_y])
    nsteps = Integer(ceil(1.5*Tend * (1/delta_x + 1/delta_y)))

    k = 50
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))

    seminormErrors = zeros(k, 1)

    for j = 1:k

        useMMS = false
        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)


        seminormErrors[j] = (SEM_Wave_2d.GradIntegral(simul, uStart - simul.uFiltered, uStart - simul.uFiltered) + SEM_Wave_2d.LpNorm(simul, uStartDer, simul.uDerFiltered, 2)^2)^(0.5)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered 

        println(string(j) * " out of " * string(k))

    end

    plt = plot(seminormErrors, yscale=:log10, label = "errors")
    savefig(plt, "seminormErrors")

end



function WedgeMain()



    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0


    #a = 1/sqrt(2)
    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    
    spaceTol = 1e-6
    systemTol = 1e-6


    #set up range of omegas
    
    
    omega = 5.0
    
    N = 8   # degree of interpolation polynomials


    c_square(x, y) = 2100.0^2

    L = maximum([xr-xl, yr-yl]) # max side length of the domain
    lambda = 2*pi*2100/omega
    N_lambda = L / lambda

    println("wavelength: " * string(lambda))
    println("length of domain in wavelengths: " * string(N_lambda))

    

    PPW = pi*(N_lambda/spaceTol)^(1/N)

    
    #lambda_min = c_min/(2*pi*omega)
    #N_lambda = L/lambda_min # upper bound on length of domain in wavelength
    #PPW = pi*(N_lambda/spaceTol)^(1/N) # upper bound on points per wavelength from rule-of-thumb


    println("(upper bound on) length of domain in wavelengths: " * string(N_lambda))
    println("Points per (shortest) wavelength: " * string(PPW))

    dofs_x = ceil(PPW * N_lambda)       # how many points do we want per dimension? points = point/wavelength * wavelength

    Kx = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!
    Ky = Integer(ceil(Kx * (5/3))) # number of elements in y-direction

    Kx = 2*Kx
    Ky = 2*Ky

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
    fVals = 100 * exp.(-5*(simul.y./1000).^2) * exp.(-5*((simul.x .- 200.0)./600).^2)'

    
    g = zeros(length(simul.y), length(simul.x))        

    println("WHI || || omega: " * string(omega) * " || " * "number of elements: " * string(Ky) * " by " * string(Kx))
    sol1, nIter = SEM_Wave_2d.WaveholtzConvHistory(simul, omega, fVals, bc, g, 200, true)



end