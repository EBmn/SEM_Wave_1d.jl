


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

    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    N = 4
    
    tol = 1e-6
    omega = 15.0

    sinc(x) = sin(2*pi*x)/(2*pi*x)
    cosc(x) = (cos(2*pi*x)-1)/(2*pi*x)
    beta(x) = sinc((x + omega)/omega) + sinc((x - omega)/omega) - 0.5*sinc(x/omega)
    gamma(x) = cosc((omega - x)/omega) - cosc((omega + x)/omega) + 0.5*cosc(x/omega)
    mu(x) = complex(beta(x), gamma(x))
    theta(x) = angle(mu(x))
    r(x) = abs(mu(x))


    c_square(x, y) = 1
    c_min = 1.0 ### has to be set properly! ###

    L = maximum([xr-xl, yr-yl]) # length of the domain

    # rule of thumb for 4th order, rounded up a bit
    dofs_x = Int(ceil(2* 1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
    numberOfElements = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!


    Kx = numberOfElements # number of elements in x-direction
    Ky = numberOfElements # number of elements in y-direction

    println("Number of elements in both directions: " * string(Kx))

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    f1 = 5
    f2 = 100
    f3 = 0.1
    fVals = f1*ones(length(simul.y), 1) * cos.(3*pi*simul.x)' + f2 * cos.(2*pi*simul.y) * ones(1, length(simul.x)) + f3 * cos.(3*pi*simul.y) * cos.(2*pi*simul.x)'
          #
    
    # the analytical solution is then
    u_sol = f1 * ones(length(simul.y), 1) * cos.(3*pi*simul.x)' / (omega^2 - 9*pi^2) + 
            f2 * cos.(2*pi*simul.y) * ones(1, length(simul.x)) / (omega^2 - 4*pi^2) + 
            f3 * cos.(3*pi*simul.y) * cos.(2*pi*simul.x)' / (omega^2 - 13*pi^2)

    println("rates: ")
    r1 = r(3*pi)
    r2 = r(2*pi)
    r3 = r(sqrt(13)*pi)
    println(r1)
    println(r2)
    println(r3)

    println("amplitudes: ")
    println(f1/(omega^2 - 9*pi^2))
    println(f2/(omega^2 - 4*pi^2))
    println(f3/(omega^2 - 13*pi^2))


    u_sol_seminorm = (SEM_Wave_2d.GradIntegral(simul, u_sol, u_sol) + SEM_Wave_2d.LpNorm(simul, u_sol, zeros(length(simul.y), length(simul.x)), 2)^2)^(0.5)


    g = zeros(length(simul.y), length(simul.x))

    ### from SEM_Wave_2d.Waveholtz:
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    spaceStep = minimum([delta_x delta_y])
    nsteps = Integer(ceil(1.5*Tend * (1/delta_x + 1/delta_y)))


    # do a lil' waveholtzin' and measure errors

    k = 80
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    realL2errors = zeros(k, 1)
    imagL2errors = zeros(k, 1)
    realH1errors = zeros(k, 1)
    imagH1errors = zeros(k, 1)
    sumErrors = zeros(k, 1)
    seminormErrors = zeros(k, 1)

    referenceLine1 = zeros(k, 1)
    referenceLine2 = zeros(k, 1)
    referenceLine3 = zeros(k, 1)

    anim = Animation()
    
    maxVal = maximum(abs.(u_sol))

    surface(simul.x, simul.y[end:-1:1], uStart, legend=:false, zlims = [-maxVal, maxVal])
    frame(anim)



    for j = 1:k

        # take a step
        useMMS = false
        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)

        # calculate errors
        #realL2errors[j] = SEM_Wave_2d.LpNorm(simul, simul.uFiltered, u_sol, 2)
        #imagL2errors[j] = SEM_Wave_2d.LpNorm(simul, simul.uDerFiltered, zeros(length(simul.y), length(simul.x)), 2)

        realH1errors[j] = SEM_Wave_2d.SobolevNorm(simul, simul.uFiltered - u_sol)
        imagH1errors[j] = SEM_Wave_2d.SobolevNorm(simul, simul.uDerFiltered)
        sumErrors[j] = (SEM_Wave_2d.SobolevNorm(simul, simul.uFiltered-u_sol)^2 + SEM_Wave_2d.SobolevNorm(simul, omega*simul.uDerFiltered)^2)^(1/2)
        seminormErrors[j] = SEM_Wave_2d.SeminormWH(simul, simul.uFiltered - u_sol, simul.uDerFiltered)


        if j == 1
            referenceLine1[1] = seminormErrors[1]
            referenceLine2[1] = seminormErrors[1]
            referenceLine3[1] = seminormErrors[1]
        else
            referenceLine1[j] = referenceLine1[j-1]*r1
            referenceLine2[j] = referenceLine2[j-1]*r2
            referenceLine3[j] = referenceLine3[j-1]*r3
        end

        # update initial conditions
        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

        # plot for fun
        surface(simul.x, simul.y[end:-1:1], simul.uFiltered, legend=:false, zlims = [-maxVal, maxVal])
        
        frame(anim)

        println(string(j) * " out of " * string(k))

    end

    println("final seminorm error: " * string(seminormErrors[end]))

    gif(anim, "waveholtzNormTest.gif", fps=10)


    
    plt = surface(simul.x, simul.y[end:-1:1], fVals, legend=:false)
    savefig(plt, "u^0_n.png")

    #=
    plt = surface(simul.x, simul.y[end:-1:1], simul.uDerFiltered, legend=:false)    
    savefig(plt, "u^1_n.png")
    =#

    #plt = plot(realL2errors, yscale=:log10, label = "L2-errors in u_0", legend=:bottomleft)
    #plot!(imagL2errors, yscale=:log10, label = "L2-errors in u_1")
    plt = plot(realH1errors, yscale=:log10, label = "H1-errors in Re(u^n)", legend=:topright)
    plot!(imagH1errors, yscale=:log10, label = "H1-errors in Im(u^n)")
    plot!(seminormErrors, yscale=:log10, label = "error in the seminorm")
    plot!(sumErrors, yscale=:log10, label = "sumErrors")
    #plot!(referenceLine1, yscale=:log10, label = "first mode rate")
    #plot!(referenceLine2, yscale=:log10, label = "second mode rate")
    #plot!(referenceLine3, yscale=:log10, label = "third mode rate")
    
    savefig(plt, "normTestErrors")

end


function main2()


        ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    N = 4
    
    tol = 1e-4
    omega = 15.0

    sinc(x) = sin(2*pi*x)/(2*pi*x)
    cosc(x) = (cos(2*pi*x)-1)/(2*pi*x)
    beta(x) = sinc((x + omega)/omega) + sinc((x - omega)/omega) - 0.5*sinc(x/omega)
    gamma(x) = cosc((omega - x)/omega) - cosc((omega + x)/omega) + 0.5*cosc(x/omega)
    mu(x) = complex(beta(x), gamma(x))
    theta(x) = angle(mu(x))
    r(x) = abs(mu(x))


    c_square(x, y) = 1
    c_min = 1.0 ### has to be set properly! ###

    L = maximum([xr-xl, yr-yl]) # length of the domain

    # rule of thumb for 4th order, rounded up a bit
    dofs_x = Int(ceil(2* 1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
    numberOfElements = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!
    #numberOfElements = 40


    Kx = numberOfElements # number of elements in x-direction
    Ky = numberOfElements # number of elements in y-direction

    println("Number of elements in both directions: " * string(Kx))

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    f1 = 5
    f2 = 100
    f3 = 0.1
    fVals = f1*ones(length(simul.y), 1) * cos.(3*pi*simul.x)' + f2 * cos.(2*pi*simul.y) * ones(1, length(simul.x)) + f3 * cos.(3*pi*simul.y) * cos.(2*pi*simul.x)'
          #
    
    # the analytical solution is then
    u_sol = f1 * ones(length(simul.y), 1) * cos.(3*pi*simul.x)' / (omega^2 - 9*pi^2) + 
            f2 * cos.(2*pi*simul.y) * ones(1, length(simul.x)) / (omega^2 - 4*pi^2) + 
            f3 * cos.(3*pi*simul.y) * cos.(2*pi*simul.x)' / (omega^2 - 13*pi^2)

    println("rates: ")
    r1 = r(3*pi)
    r2 = r(2*pi)
    r3 = r(sqrt(13)*pi)
    println(r1)
    println(r2)
    println(r3)

    println("amplitudes: ")
    println(f1/(omega^2 - 9*pi^2))
    println(f2/(omega^2 - 4*pi^2))
    println(f3/(omega^2 - 13*pi^2))


    u_sol_seminorm = (SEM_Wave_2d.GradIntegral(simul, u_sol, u_sol) + SEM_Wave_2d.LpNorm(simul, u_sol, zeros(length(simul.y), length(simul.x)), 2)^2)^(0.5)


    g = zeros(length(simul.y), length(simul.x))

    ### from SEM_Wave_2d.Waveholtz:
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    spaceStep = minimum([delta_x delta_y])
    nsteps = Integer(ceil(1.5*Tend * (1/delta_x + 1/delta_y)))


    # do a lil' waveholtzin' and measure errors

    k = 80
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    realL2errors = zeros(k, 1)
    imagL2errors = zeros(k, 1)
    realH1errors = zeros(k, 1)
    imagH1errors = zeros(k, 1)
    sumErrors = zeros(k, 1)
    seminormErrors = zeros(k, 1)

    referenceLine1 = zeros(k, 1)
    referenceLine2 = zeros(k, 1)
    referenceLine3 = zeros(k, 1)

    anim = Animation()
    
    maxVal = maximum(abs.(u_sol))

    surface(simul.x, simul.y[end:-1:1], uStart, legend=:false, zlims = [-maxVal, maxVal])
    frame(anim)



    for j = 1:k

        # take a step
        useMMS = false
        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)

        # calculate errors
        #realL2errors[j] = SEM_Wave_2d.LpNorm(simul, simul.uFiltered, uStart, 2)
        
        #imagL2errors[j] = SEM_Wave_2d.LpNorm(simul, simul.uDerFiltered, uStartDer, 2)

        realH1errors[j] = SEM_Wave_2d.SobolevNorm(simul, simul.uFiltered-uStart)
        imagH1errors[j] = SEM_Wave_2d.SobolevNorm(simul, simul.uDerFiltered - uStartDer)
        sumErrors[j] = (SEM_Wave_2d.SobolevNorm(simul, simul.uFiltered-uStart)^2 + SEM_Wave_2d.SobolevNorm(simul, omega*(simul.uDerFiltered - uStartDer))^2)^(1/2)
        seminormErrors[j] = SEM_Wave_2d.SeminormWH(simul, simul.uFiltered - uStart, simul.uDerFiltered-uStartDer)

        if j == 1
            referenceLine1[1] = seminormErrors[1]
            referenceLine2[1] = seminormErrors[1]
            referenceLine3[1] = seminormErrors[1]
        else
            referenceLine1[j] = referenceLine1[j-1]*r1
            referenceLine2[j] = referenceLine2[j-1]*r2
            referenceLine3[j] = referenceLine3[j-1]*r3
        end

        # update initial conditions
        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

        # plot for fun
        surface(simul.x, simul.y[end:-1:1], simul.uFiltered, legend=:false, zlims = [-maxVal, maxVal])
        
        frame(anim)

        println(string(j) * " out of " * string(k))

    end

    println("final seminorm error: " * string(seminormErrors[end]))

    gif(anim, "waveholtzNormTest.gif", fps=10)


    
    plt = surface(simul.x, simul.y[end:-1:1], fVals, legend=:false)
    savefig(plt, "u^0_n.png")

    #=
    plt = surface(simul.x, simul.y[end:-1:1], simul.uDerFiltered, legend=:false)    
    savefig(plt, "u^1_n.png")
    =#

    #plt = plot(realL2errors, yscale=:log10, label = "L2-errors in u_0", legend=:bottomleft)
    #plot!(imagL2errors, yscale=:log10, label = "L2-errors in u_1")
    plt = plot(realH1errors, yscale=:log10, label = "H1-errors in Re(u^n)", legend=:bottomleft)
    plot!(imagH1errors, yscale=:log10, label = "H1-errors in Im(u^n)")
    plot!(seminormErrors, yscale=:log10, label = "error in the seminorm")
    plot!(sumErrors, yscale=:log10, label = "sumErrors")
    #plot!(referenceLine1, yscale=:log10, label = "first mode rate")
    #plot!(referenceLine2, yscale=:log10, label = "second mode rate")
    #plot!(referenceLine3, yscale=:log10, label = "third mode rate")
    
    savefig(plt, "normTestDiffErrors")




end