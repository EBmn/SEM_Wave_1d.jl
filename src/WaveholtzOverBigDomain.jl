using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Base.Threads
using ThreadTools

include("SEM_Wave_1d_cleaned.jl")
using .SEM_Wave_1d

using LaTeXStrings
using MAT

include("MMS.jl")
using .MMS


function main()



    xl = -100.0
    xr = 100.0

    spaceTol = 1e-1
    g = [0.0, 0.0]
    omega = 3.13

    lambda = 1/(omega/(2*pi))
    N_lambda = (xr-xl)/lambda
    
    N = 6                                                   # order of method

    PPW = pi*(N_lambda / spaceTol)^(1/N)                    # points per wavelength
    dofs_x = ceil(PPW * N_lambda)                           # how many points do we want in space?

    K = Int(ceil((dofs_x - 1)/N))                           # the corresponding number of elements!
    
    
    println("Using $K elements, corresponding to about $PPW points per wavelengths")

    

    bc = [1/sqrt(2), sqrt(1/2)]                             # impedance conditions

    c_square(x) = 1.0                                       # constant wavespeed



    simul = SEM_Wave_1d.SEM_Wave([xl, xr], K, N, c_square)


    simul.omega = omega


    simul.bc = bc
    
    
    fVals = exp.(-10*((simul.x)).^2)


    # find Helmholtz solution using the matrices:

    # the wave equation is MU'' + BU' + AU = F, corresponding to a Helmholtz equation of (-\omega^2M + i\omegaB + A)U = F
    # in the terminology of the SEM code we have M = simul.M, B = (alpha/beta)*simul.M_b, and L*U = LaplaceTerm*U / ((simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)))

    M = spdiagm(simul.M)
    B = spdiagm((simul.bc[1]/simul.bc[2]) * simul.M_b)

    # we will now reconstruct the matrix A from the action of LaplaceTerm
    A = zeros(length(simul.x), length(simul.x))

    println("$(length(simul.x)) points in space")


    for j = 1:length(simul.x)
        
        println(j)

        e_j = zeros(length(simul.x))
        e_j[j] = 1.0
        

        Ae_j = (simul.M + 0.5*(simul.bc[1]/simul.bc[2])*simul.timestep*simul.M_b) .* (SEM_Wave_1d.LaplaceTerm(simul, e_j) ./ simul.timestep^2)
        #Ae_j = SEM_Wave_1d.LaplaceTerm(simul, e_j) / simul.timestep^2

        #println(length(SEM_Wave_1d.LaplaceTerm(simul, e_j)))
        #println(length(Ae_j))

        A[:, j] = Ae_j

    end

    A = sparse(A)
    #plt = surface(simul.x, simul.x, A)
    #plt = heatmap(log10.(abs.(A)))
    #plt = plot(simul.x[2:end-1], (A*(sin.(simul.x)))[2:end-1])
    #savefig(plt, "AmatrixTest")
    



    HelmOp = (-simul.omega^2 * M + (1im*simul.omega)*B + A)

    HelmOp = (simul.omega^2 * M + (1im*simul.omega)*B + A)


    Btest = (1im*simul.omega)*B
    
    println("simul.omega is $(simul.omega)")

    println("****")
    #println(Matrix(Btest))
    println(maximum(abs.(Btest)))
    println("****")
    #println("The Helmholtz operator has condition number $(cond(HelmOp))")

    F = M*fVals;

    U_helm = HelmOp\F;

    #oFs = length(simul.x)
    #U_direct, history = gmres(HelmOp, fVals, log=true, verbose=true, reltol=1e-6,  restart=DoFs)


    
    plt = plot(simul.x, real(U_helm))
    savefig(plt, "DirectSolveTest")
    #plot(simul.x, fVals)
    
    #=
    U_whi, nIter = SEM_Wave_1d.Waveholtz(simul, fVals, omega, bc, g, 1e-0)

    println("WaveHoltz converged in $nIter iterations")
    
    plt = plot(simul.x, U_whi)
    savefig("BigDomainWhi")
    =#

    




    #Waveholtz-specific parameters for the wave solver
    simul.fVals = -fVals
    simul.omega = omega
    simul.Tend = 2*pi/omega
    

    #appropriate parameters for the wave solver:
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
    simul.timestep = simul.Tend/simul.nsteps

    #starting guess
    uStart = zeros(length(simul.x))
    uStartDer = zeros(length(simul.x))

    
    animate = false     #we do not want any animations of the wave eq solutions


    simul.useWaveholtz = true


    nIter = 100

    errors = zeros(nIter)

    bracket(x) = (1 + x^2)^(-1.6)

    bracketVec = bracket.(simul.x)

    for j = 1:nIter
        
        println("iteration: $j")
        SEM_Wave_1d.Simulate(simul, uStart, uStartDer, simul.Tend, simul.nsteps, simul.fVals, simul.omega, simul.bc, simul.g)

        uStart = simul.uFiltered
        #uStartDer = simul.uDerFiltered
        uStartDer = zeros(length(simul.x))

        #println(norm(uStart))
        #println(norm((real(U_helm)-uStart) .* bracketVec))
        println(SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, uStart .* bracketVec, 2))

        #errors[j] = SEM_Wave_1d.LpNorm(simul, real(U_helm), uStart, 2)
        errors[j] = SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, uStart .* bracketVec, 2)
        
        #errors[j] = norm(real(U_helm)-uStart)

    end

    println(errors)
    #plt = plot(LinRange(1, nIter, nIter), errors .+ 1 ./ (sqrt.(LinRange(1, nIter, nIter))), xscale=:log10, yscale=:log10)
    plt = plot(LinRange(1, nIter, nIter), errors, xscale=:log10, yscale=:log10)
    plot!(LinRange(1, nIter, nIter), errors[1] ./ (sqrt.(LinRange(1, nIter, nIter))), xscale=:log10, yscale=:log10, label=L"O(n^{-1/2})")

    savefig(plt, "tester")

    plt = plot(simul.x, uStart)
    savefig(plt, "Waveholtz_Iterate.pdf")


    #=
    matwrite("data.mat", Dict(
    "M" => M,
    "A" => A,
    "B" => B,
    "HelmOp" => HelmOp,
    "U_dir" => U_direct,
    "F_julia" => F,
    "U_whi" => U_whi
    ))
    =#


end


function differentWeightsPlot()


    nIter = 500                                             # this many iterations means the support of the iterates is about [-nIter * T. nIter*T]

    println(nIter)

    #spaceTol = 1e-8 seems not to be worth it!
    spaceTol = 1e-4 # seems to be good enough?
    g = [0.0, 0.0]
    omega = 8.13
    
    #omega = 2*omega

    xl = -50.0
    xr = 50.0



    lambda = 1/(omega/(2*pi))                               # wavelength
    N_lambda = (xr-xl)/lambda                               # length of domain in wavelengths
    
    N = 8                                                   # order of method

    PPW = pi*(N_lambda / spaceTol)^(1/N)                    # points per wavelength
    dofs_x = ceil(PPW * N_lambda)                           # how many points do we want in space?

    K = Int(ceil((dofs_x - 1)/N))                           # the corresponding number of elements!


    #K = Int(2*K)


    println("Using $K elements, corresponding to about $PPW points per wavelengths")



    bc = [1/sqrt(2), sqrt(1/2)]                             # impedance conditions

    #c_square(x) = 1.0                                       # constant wavespeed
    c_square(x) = 1.0 - 0.9*exp(-5*(x-1.0)^2) - 0.9*exp(-5*(x+1.0)^2)  # nonconstant wavespeed


    simul = SEM_Wave_1d.SEM_Wave([xl, xr], K, N, c_square)

    simul.omega = omega
    simul.bc = bc


    fVals = exp.(-10*((simul.x)).^2)                        # forcing


    # find Helmholtz solution using the matrices:

    # the wave equation is MU'' + BU' + AU = F, corresponding to a Helmholtz equation of (-\omega^2M + i\omegaB + A)U = F
    # in the terminology of the SEM code we have M = simul.M, B = (alpha/beta)*simul.M_b, and L*U = LaplaceTerm*U / ((simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)))

    M = spdiagm(simul.M)
    B = spdiagm((simul.bc[1]/simul.bc[2]) * simul.M_b)

    # we will now reconstruct the matrix A from the action of LaplaceTerm
    A = zeros(length(simul.x), length(simul.x))

    println("$(length(simul.x)) points in space")


    for j = 1:length(simul.x)
        
        println("generating row $j of A:")

        e_j = zeros(length(simul.x))
        e_j[j] = 1.0
        
        Ae_j = (simul.M + 0.5*(simul.bc[1]/simul.bc[2])*simul.timestep*simul.M_b) .* (SEM_Wave_1d.LaplaceTerm(simul, e_j) ./ simul.timestep^2)

        A[:, j] = Ae_j

    end

    A = sparse(A)
    

    H = (simul.omega^2 * M + (1im*simul.omega)*B + A)       # construct the Helmholtz operator


    Btest = (1im*simul.omega)*B
    

    F = M*fVals;                                            # the matrix F
    U_helm = H\F;                                           # Direct solution


    plt = plot(simul.x, real(U_helm), label="real part")
    plot!(simul.x, imag(U_helm), label="imaginary part")
    savefig(plt, "DirectSolveTest.png")

    ############### WaveHoltz iterations ###############


    #Waveholtz-specific parameters for the wave solver
    simul.fVals = -fVals
    simul.omega = omega
    simul.Tend = 2*pi/omega
    

    #appropriate parameters for the wave solver:
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
    simul.timestep = simul.Tend/simul.nsteps

    #starting guess
    uStart = zeros(length(simul.x))
    uStartDer = zeros(length(simul.x))

    
    animate = false                                         # we do not want any animations of the wave eq solutions

    simul.useWaveholtz = true                               # we do want to use Waveholtz 
    
    
    sVals = [0.5, 1, 1.5, 2, 3, 4, 5, 6]

    errors = zeros(nIter, length(sVals))
    uNorms = zeros(nIter, length(sVals))



    #plt = plot(LinRange(1, nIter, nIter), 1 ./ (sqrt.(LinRange(1, nIter, nIter))), xscale=:log10, yscale=:log10, label=L"O(n^{-1/2})")
    plt = plot()

    

    # reset starting guess and filtered solution
    uStart = zeros(length(simul.x))
    uStartDer = zeros(length(simul.x))
    simul.uFiltered = zeros(length(simul.x))
    


    for j = 1:nIter
        
        println("iteration: $j")
        SEM_Wave_1d.Simulate(simul, uStart, uStartDer, simul.Tend, simul.nsteps, simul.fVals, simul.omega, simul.bc, simul.g)

        uStart = simul.uFiltered
        #uStartDer = simul.uDerFiltered
        uStartDer = zeros(length(simul.x))

        for k = 1:length(sVals)

            s = sVals[k]
            bracket(x) = (1 + x^2)^(-s/2)
            bracketVec = bracket.(simul.x)
            errors[j, k] = SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, uStart .* bracketVec, 2)
            uNorms[j, k] = SEM_Wave_1d.LpNorm(simul, uStart  .* bracketVec, zeros(length(simul.x)), 2)

        end

        
        
    end


    for k = 1:length(sVals)
        
        s = sVals[k]
        legendString = L"s = " * string(s)
        plot!(LinRange(1, nIter, nIter), errors[:, k], xscale=:log10, yscale=:log10, label=legendString)
        #plot!(LinRange(1, nIter, nIter), errors[:, k] .+ 1 ./ (sqrt.(LinRange(1, nIter, nIter))), xscale=:log10, yscale=:log10, label=legendString)

    end

    plot!(LinRange(1, nIter, nIter), errors[1, 1] ./ (sqrt.(LinRange(1, nIter, nIter))), xscale=:log10, yscale=:log10, label=L"O(n^{-1/2})", linestyle=:dash)

    savefig(plt, "PlotDifferentS.pdf")
    


    plt = plot(simul.x, simul.uFiltered)
    #plt = plot(simul.x, (simul.uFiltered - real(U_helm)) ./ real(U_helm))
    savefig(plt, "FinalIterate.pdf")



    matwrite("data.mat", Dict(
    "M" => M,
    "A" => A,
    "B" => B,
    "H" => H,
    "U_dir" => U_helm,
    "F" => F,
    "U_whi" => simul.uNow,
    "uNorms" => uNorms
    ))



end




function differentOmegasPlot()


    #waveholtzTol = 1
    waveholtzTol = 5e-3
    spaceTol = 1e-2
    
    g = [0.0, 0.0]

    xl = -300.0                                             # size of domain to avoid reacing the boundary with WaveHoltz 
    xr = 300.0

    omegaMax = 15
    omegas = LinRange(2, omegaMax, omegaMax-1)

    lambda = 1/(omegaMax/(2*pi))                               # wavelength
    N_lambda = (xr-xl)/lambda                               # length of domain in wavelengths
    
    N = 8                                                   # order of method

    PPW = pi*(N_lambda / spaceTol)^(1/N)                    # points per wavelength
    dofs_x = ceil(PPW * N_lambda)                           # how many points do we want in space?

    K = Int(ceil((dofs_x - 1)/N))                           # the corresponding number of elements!


    nIters = zeros(length(omegas))
    

    sVals = [0.5, 1, 1.5, 2, 3, 4, 5, 6]                    # the different s parameters to use for norms

    bc = [1/sqrt(2), sqrt(1/2)]                             # impedance conditions

    c_square(x) = 1.0                                       # constant wavespeed

    simul = SEM_Wave_1d.SEM_Wave([xl, xr], K, N, c_square)
    simul.bc = bc

    println("Using $K elements, i.e. $(length(simul.x)) points in space")


    ###################### construct the matrices ######################
    M = spdiagm(simul.M)
    B = spdiagm((simul.bc[1]/simul.bc[2]) * simul.M_b)

    # we will now reconstruct the matrix A from the action of LaplaceTerm
    A = zeros(length(simul.x), length(simul.x))

    println("$(length(simul.x)) points in space")


    
    simul.omega = 1.0; # nonzero placeholder
    
    for j = 1:length(simul.x)
        
        println("generating row $j of A:")

        e_j = zeros(length(simul.x))
        e_j[j] = 1.0
        
        Ae_j = (simul.M + 0.5*(simul.bc[1]/simul.bc[2])*simul.timestep*simul.M_b) .* (SEM_Wave_1d.LaplaceTerm(simul, e_j) ./ simul.timestep^2)

        A[:, j] = Ae_j

    end


    A = sparse(A)


    maxIter = 1e4

    uSolsReal = zeros(length(simul.x), length(omegas))
    uSolsImag = zeros(length(simul.x), length(omegas))
    uSolNorms = zeros(length(omegas), Integer(maxIter), length(sVals))

    

    s = 2.0

    for k = 1:length(omegas)

        omega = omegas[k]
        simul.omega = omega

        S(x) = exp(-10*(x).^2)
        f(x) = omega^2 * S(omega * x)
        fVals = f.(simul.x)                                     # forcing
        F = M*fVals;                                            # the matrix F

        # find Helmholtz solution using the matrices:

        # the wave equation is MU'' + BU' + AU = F, corresponding to a Helmholtz equation of (-\omega^2M + i\omegaB + A)U = F
        # in the terminology of the SEM code we have M = simul.M, B = (alpha/beta)*simul.M_b, and L*U = LaplaceTerm*U / ((simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)))



        H = (simul.omega^2 * M + (1im*simul.omega)*B + A)       # construct the Helmholtz operator

        U_helm = H\F;                                           # Direct solution
        uSolsReal[:, k] = real(U_helm)
        uSolsImag[:, k] = imag(U_helm)



    


        plt = plot(simul.x, real(U_helm), label="real part")
        plot!(simul.x, imag(U_helm), label="imaginary part")
        savefig(plt, "DirectSolveTest.png")

        ############### WaveHoltz solves ###############
        #Waveholtz-specific parameters for the wave solver
        simul.fVals = -fVals
        simul.omega = omega
        simul.Tend = 2*pi/omega
        

        # appropriate parameters for the wave solver:
        delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
        simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
        simul.timestep = simul.Tend/simul.nsteps

        # starting guess
        uStart = zeros(length(simul.x))
        uStartDer = zeros(length(simul.x))

        # reset filtered solutions for WaveHoltz
        simul.uFiltered = zeros(length(simul.x))
        simul.uDerFiltered = zeros(length(simul.x))
        
        animate = false                                         # we do not want any animations of the wave eq solutions

        simul.useWaveholtz = true                               # we do want to use Waveholtz 

        
        bracket(x) = (1 + x^2)^(-s/2)
        bracketVec = bracket.(simul.x)

        res = Inf
        
        nIter = 0

        while res > waveholtzTol
            
            nIter = nIter + 1

            SEM_Wave_1d.Simulate(simul, uStart, uStartDer, simul.Tend, simul.nsteps, simul.fVals, simul.omega, simul.bc, simul.g)

            uStart = simul.uFiltered
            uStartDer = zeros(length(simul.x))

            res = SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, uStart .* bracketVec, 2)# / (SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, zeros(length(simul.x)), 2))
            
            println("iteration: $nIter || omega = $omega || residual = $res" )

            plt = plot(simul.x, uStart)
            savefig(plt, "iterateTest.png")
            
            if nIter > maxIter
                println("maximum iteration count reached!")
                res = 0.0
            end


            for j = 1:length(sVals)

                s = sVals[j]
                bracket(x) = (1 + x^2)^(-s/2)
                bracketVec = bracket.(simul.x)
                uSolNorms[k, nIter, j] = SEM_Wave_1d.LpNorm(simul, real(U_helm)  .* bracketVec, uStart .* bracketVec, 2)

            end


            #println("weighted norm of Re(U): $(SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, zeros(length(simul.x)), 2)))")

        end

        nIters[k] = nIter

    end


    plt = plot(omegas, nIters, xscale=:log10, yscale=:log10, label="iterations")
    #plot!(omegas, (nIters[end]/(omegas[end])^(2*s-1)) .* omegas.^(2*s-1), xscale=:log10, yscale=:log10, label=L"O(\omega^{2s-1})", linestyle=:dash)
    savefig(plt, "PlotDifferentOmegas.pdf")


    matwrite("differentOmegaData.mat", Dict(
    "M" => M,
    "A" => A,
    "B" => B,
    "omegas" => omegas,
    "sVals" => sVals,
    "U_sols_real" => uSolsReal,
    "U_sols_imag" => uSolsImag,
    "U_whi" => simul.uNow,
    "uNorms" => uSolNorms,
    "xVals" => simul.x
    ))





end



function differentOmegasFixedIter()


    #waveholtzTol = 1
    
    #spaceTol = 1e-4
    spaceTol = 1e-3
    
    g = [0.0, 0.0]

    xl = -50.0                                             # size of domain to avoid reacing the boundary with WaveHoltz 
    xr = 50.0


    # we have an exterior error bounded above by sup_{xl<x<xr}(U) * \int_{x>xr} \langle x\rangle^{-2s}dx,
    # with |xl| = |xr| = 100 and s = 3/2 this comes out to be approximately sup_{xl<x<xr}(U)*5e-5
    # for the forcing we have chosen below, we have a supremum of the order of 1e-1,
    # => should have exterior error on the order of 1e-6


    maxIter = 50

    omegaMax = 50
    omegas = LinRange(10, omegaMax, omegaMax)

    lambda = 1/(omegaMax/(2*pi))                               # wavelength
    N_lambda = (xr-xl)/lambda                               # length of domain in wavelengths
    
    N = 8                                                   # order of method

    PPW = pi*(N_lambda / spaceTol)^(1/N)                    # points per wavelength
    dofs_x = ceil(PPW * N_lambda)                           # how many points do we want in space?

    K = Int(ceil((dofs_x - 1)/N))                           # the corresponding number of elements!


    nIters = zeros(length(omegas))
    

    sVals = [0.5, 1, 1.5, 2, 3, 4, 5, 6]                    # the different s parameters to use for norms

    bc = [1/sqrt(2), sqrt(1/2)]                             # impedance conditions

    c_square(x) = 1.0                                       # constant wavespeed

    simul = SEM_Wave_1d.SEM_Wave([xl, xr], K, N, c_square)
    simul.bc = bc

    println("Using $K elements, i.e. $(length(simul.x)) points in space")


    ###################### construct the matrices ######################
    M = spdiagm(simul.M)
    B = spdiagm((simul.bc[1]/simul.bc[2]) * simul.M_b)

    # we will now reconstruct the matrix A from the action of LaplaceTerm
    A = zeros(length(simul.x), length(simul.x))

    println("$(length(simul.x)) points in space")


    uDirNorms = zeros(length(omegas))


    
    simul.omega = 1.0; # nonzero placeholder
    
    for j = 1:length(simul.x)
        
        println("generating row $j of A:")

        e_j = zeros(length(simul.x))
        e_j[j] = 1.0
        
        Ae_j = (simul.M + 0.5*(simul.bc[1]/simul.bc[2])*simul.timestep*simul.M_b) .* (SEM_Wave_1d.LaplaceTerm(simul, e_j) ./ simul.timestep^2)

        A[:, j] = Ae_j

    end


    A = sparse(A)


    

    uSolsReal = zeros(length(simul.x), length(omegas))
    uSolsImag = zeros(length(simul.x), length(omegas))
    uSolNorms = zeros(length(omegas), Integer(maxIter), length(sVals))

    

    s = 2.0

    for k = 1:length(omegas)

        omega = omegas[k]
        simul.omega = omega

        S(x) = exp(-10*(x).^2)
        f(x) = omega^2 * S(omega * x)
        fVals = f.(simul.x)                                     # forcing
        F = M*fVals;                                            # the matrix F

        # find Helmholtz solution using the matrices:

        # the wave equation is MU'' + BU' + AU = F, corresponding to a Helmholtz equation of (-\omega^2M + i\omegaB + A)U = F
        # in the terminology of the SEM code we have M = simul.M, B = (alpha/beta)*simul.M_b, and L*U = LaplaceTerm*U / ((simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)))



        H = (simul.omega^2 * M + (1im*simul.omega)*B + A)       # construct the Helmholtz operator

        U_helm = H\F;                                           # Direct solution
        uSolsReal[:, k] = real(U_helm)
        uSolsImag[:, k] = imag(U_helm)


        #uDirNorms = SEM_Wave_1d.LpNorm(simul, (sqrt.(real(U_helm) + imag(U_helm))) .* bracketVec, zeros(length(simul.x)), 2)

    


        plt = plot(simul.x, real(U_helm), label="real part")
        plot!(simul.x, imag(U_helm), label="imaginary part")
        savefig(plt, "DirectSolveTest.png")

        ############### WaveHoltz solves ###############
        #Waveholtz-specific parameters for the wave solver
        simul.fVals = -fVals
        simul.omega = omega
        simul.Tend = 2*pi/omega
        

        # appropriate parameters for the wave solver:
        delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
        simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
        simul.timestep = simul.Tend/simul.nsteps

        # starting guess
        uStart = zeros(length(simul.x))
        uStartDer = zeros(length(simul.x))

        # reset filtered solutions for WaveHoltz
        simul.uFiltered = zeros(length(simul.x))
        simul.uDerFiltered = zeros(length(simul.x))
        
        animate = false                                         # we do not want any animations of the wave eq solutions

        simul.useWaveholtz = true                               # we do want to use Waveholtz 

        
        bracket(x) = (1 + x^2)^(-s/2)
        bracketVec = bracket.(simul.x)

        res = Inf


        for nIter = 1:maxIter

            SEM_Wave_1d.Simulate(simul, uStart, uStartDer, simul.Tend, simul.nsteps, simul.fVals, simul.omega, simul.bc, simul.g)

            uStart = simul.uFiltered
            uStartDer = zeros(length(simul.x))

            res = SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, uStart .* bracketVec, 2)# / (SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, zeros(length(simul.x)), 2))
            
            println("iteration: $nIter || omega = $omega || residual = $res" )

            plt = plot(simul.x, uStart)
            savefig(plt, "iterateTest.png")
            
            if nIter > maxIter
                println("maximum iteration count reached!")
                res = 0.0
            end


            for j = 1:length(sVals)

                s = sVals[j]
                bracket(x) = (1 + x^2)^(-s/2)
                bracketVec = bracket.(simul.x)
                uSolNorms[k, nIter, j] = SEM_Wave_1d.LpNorm(simul, real(U_helm)  .* bracketVec, uStart .* bracketVec, 2)

            end


            #println("weighted norm of Re(U): $(SEM_Wave_1d.LpNorm(simul, real(U_helm) .* bracketVec, zeros(length(simul.x)), 2)))")

        end

    end


    #plt = plot(omegas, uSolNorms[], xscale=:log10, yscale=:log10, label="iterations")
    #plot!(omegas, (nIters[end]/(omegas[end])^(2*s-1)) .* omegas.^(2*s-1), xscale=:log10, yscale=:log10, label=L"O(\omega^{2s-1})", linestyle=:dash)
    #savefig(plt, "PlotFixedIter.pdf")


    matwrite("differentOmegaFixedIter.mat", Dict(
    "M" => M,
    "A" => A,
    "B" => B,
    "omegas" => omegas,
    "sVals" => sVals,
    "U_sols_real" => uSolsReal,
    "U_sols_imag" => uSolsImag,
    "U_whi" => simul.uNow,
    "uNorms" => uSolNorms,
    "xVals" => simul.x
    ))





end





function WaveTest()



    xl = -10.0
    xr = 10.0

    g = [0.0, 0.0]
    omega = 8.13

    minDx = 1e-1                                            # tolerance on element size
    K = Integer(ceil((xr-xl)/minDx))
    N = 4                                                   # order of method

    bc = [1/sqrt(2), sqrt(1/2)]                             # impedance conditions

    c_square(x) = 1.0                                       # constant wavespeed

    simul = SEM_Wave_1d.SEM_Wave([xl, xr], K, N, c_square)
    
    fVals = 20*exp.(-((simul.x.-0.5)/0.1).^2)

    uStart = 5*exp.(-((simul.x.-0.5)/0.1).^2)
    uStartDer = zeros(length(simul.x))
    


    nsteps = 1500
    Tend = 15.0
    
    SEM_Wave_1d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, true)




end
