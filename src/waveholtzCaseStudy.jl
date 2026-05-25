


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
    
    tol = 1e-2
    #omega = 9.225806451612902
    omega = 8.0

    c_square(x, y) = 1
    c_min = 1.0 ### has to be set properly! ###

    L = maximum([xr-xl, yr-yl]) # length of the domain

    # rule of thumb for 4th order, rounded up a bit
    dofs_x = Int(ceil(4*1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
    numberOfElements = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!


    Kx = numberOfElements # number of elements in x-direction
    Ky = numberOfElements # number of elements in y-direction

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    fVals = ones(length(simul.y), 1) * cos.(3*pi*simul.x)' + 150 * ones(length(simul.y), 1) * cos.(pi*simul.x)'
    # this means that the solution is precisely cos.(omega*x)/(omega^2 - 9*pi^2)
    
    u_sol = ones(length(simul.y), 1) * cos.(3*pi*simul.x)' / (omega^2 - 9*pi^2)
    u_sol_seminorm = (SEM_Wave_2d.GradIntegral(simul, u_sol, u_sol) + SEM_Wave_2d.LpNorm(simul, u_sol, zeros(length(simul.y), length(simul.x)), 2)^2)^(0.5)
    

    #= 
    this means that the solution to the wave equation with zero initial condition is precisely 
    (1/(omega^2 - 9*pi^2)) * (cos(omega*t) - cos(3*pi*t)) * cos(3*pi*x) + cos(4*pi*t)*cos(4*pi*x) + (1/(5*pi)) * sin(5*pi*t)*cos(5*pi*x)
    =#


    g = zeros(length(simul.y), length(simul.x))        

    ### from SEM_Wave_2d.Waveholtz:
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    spaceStep = minimum([delta_x delta_y])
    nsteps = Integer(ceil(1.5*Tend * (1/delta_x + 1/delta_y)))



    # do a lil' waveholtzin' and measure errors

    #u_0, u_1, history = SEM_Wave_2d.WaveholtzGMRESnew(simul, omega, fVals, bc, g, 1e-9)
    


    k = 50
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    errors = zeros(k, 1)
    errors2 = zeros(k, 1)
    errors3 = zeros(k, 1)
    impliedBetas = zeros(k-1, 1)
    impliedGammas = zeros(k-1, 1)
    
    predictedRealErrors = zeros(k, 1)
    predictedImagErrors = zeros(k, 1)
    referenceLine = zeros(k, 1)
    referenceLine2 = zeros(k, 1)

    sinc(x) = sin(2*pi*x)/(2*pi*x)
    cosc(x) = (cos(2*pi*x)-1)/(2*pi*x)
    beta = sinc((3*pi + omega)/omega) + sinc((3*pi - omega)/omega) - 0.5*sinc(3*pi/omega)
    gamma = cosc((omega - 3*pi)/omega) - cosc((omega +3*pi)/omega) + 0.5*cosc(3*pi/omega)
    mu = complex(beta, gamma)
    theta = angle(mu)
    r = abs(mu)

    anim = Animation()
    surface(simul.x, simul.y[end:-1:1], uStart, legend=:false, zlims = [-1, 1])
    frame(anim)


    for j = 1:k

        # found by hand:
        predictedRealErrors[j] = abs((1/(9*pi^2 - omega^2)) * 0.25 * r^j * cos(j*theta))
        predictedImagErrors[j] = abs((3*pi/(9*pi^2 - omega^2)) * 0.25 * r^j * sin(j*theta))
        referenceLine[j] = r^j * u_sol_seminorm

        useMMS = false
        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)

        errors[j] = SEM_Wave_2d.LpNorm(simul, simul.uFiltered, uStart, 2)
        errors2[j] = SEM_Wave_2d.LpNorm(simul, simul.uDerFiltered, uStartDer, 2)
        errors3[j] = (SEM_Wave_2d.GradIntegral(simul, simul.uFiltered - uStart, simul.uFiltered - uStart) + SEM_Wave_2d.LpNorm(simul, simul.uDerFiltered, uStartDer, 2)^2)^(0.5)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered # uStartDer should approach zero...

        #uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))

        #surface(simul.x, simul.y[end:-1:1], u_sol - simul.uFiltered, legend=:false, zlims = [-1, 1])
        surface(simul.x, simul.y[end:-1:1], simul.uFiltered, legend=:false, zlims = [-1, 1])
        #surface(simul.x, simul.y[end:-1:1], simul.uDerFiltered, legend=:false, zlims = [-3, 3])
        frame(anim)

        #errors[j] = SEM_Wave_2d.LpNorm(simul, simul.uFiltered, u_sol, 2)
        #errors2[j] = SEM_Wave_2d.LpNorm(simul, simul.uDerFiltered, zeros(length(simul.y), length(simul.x)), 2)

        #errors3[j] = (SEM_Wave_2d.GradIntegral(simul, simul.uFiltered - u_sol, simul.uFiltered - u_sol) + SEM_Wave_2d.LpNorm(simul, simul.uDerFiltered, zeros(length(simul.y), length(simul.x)), 2)^2)^(0.5)
        

        if j > 1
            impliedBetas[j-1] = errors[j]/errors[j-1] 
            impliedGammas[j-1] = errors2[j]/errors2[j-1] 
        end

        println(string(j) * " out of " * string(k))

    end

    println("final error: " * string(errors3[end]))

    gif(anim, "waveholtzCaseStudy.gif", fps=10)


    plt = surface(simul.x, simul.y[end:-1:1], simul.uFiltered, legend=:false)
    #plt = heatmap(simul.x, simul.y, simul.uFiltered, legend=:false)
    savefig(plt, "u^0_n.png")

    plt = surface(simul.x, simul.y[end:-1:1], simul.uDerFiltered, legend=:false)    
    savefig(plt, "u^1_n.png")


    println("r_k: " * string(r))
    println("theta_k: " * string(theta))



    plt = plot(errors, yscale=:log10, label = "errors in u_0")
    #plot!(referenceLine, yscale=:log10)
    #plot!(predictedRealErrors, yscale=:log10, label = "ex L^2 errors in u_0")
    plot!(errors2, yscale=:log10, label = "errors in u_1")
    #plot!(predictedImagErrors, yscale=:log10, label = "ex L^2 errors in u_1")
    plot!(errors3, yscale=:log10, label = "error in the ||| . |||-norm", legend=:bottomleft)
    savefig(plt, "caseStudyErrors")

    plt = plot(impliedBetas, label = "numerical beta")
    plot!(ones(k-1, 1)*0.997054653246982, label = "analytical beta")
    plot!(impliedGammas, label = "numerical gamma")
    plot!(ones(k-1, 1)*0.067658235513621, label = "analytical gamma")

    savefig(plt, "impliedBetas")


end