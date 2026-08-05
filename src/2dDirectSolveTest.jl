
using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using IterativeSolvers
using TimerOutputs

include("SEM_Wave_2d_Updated.jl")
using .SEM_Wave_2d_Updated
using LaTeXStrings
using MAT

include("MMS.jl")
using .MMS


function main()
# domain
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0

    omega = 5.1


    N = 8   # degree of interpolation polynomials
    Kx = 4 # number of elements in x-direction
    Ky = 4 # number of elements in y-direction

    heaviside(x) = 0.5 * (sign(x) + 1)
    box(x, a, b) = heaviside(x-a) - heaviside(x-b)

    c_square(x, y) = 1# - 0.9*box(x, 0.2, 0.5)*box(y, 0.3, 0.6) - 0.9*box(x, -0.5, -0.2)*box(y, -0.6, -0.3)

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    nx = length(simul.x)
    ny = length(simul.y)
    
    f(x, y) = exp.(-((y .- 0.2*(yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
    fVals = f.(simul.x, simul.y')

    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] #impedance all around

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    

    

    a = simul.stepCoeffs[1]
    b = simul.stepCoeffs[2]
    
    n = Int(length(simul.x) * length(simul.y))

    H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)


    uSolVec = H\F # we could also solve this system using e.g. gmres, since we have the matrices!
    
    uSol = reshape(uSolVec, ny, nx)

    plt = heatmap(simul.y, simul.x, real(uSol))

    savefig(plt, "directSolReal.pdf")

    plt = heatmap(simul.y, simul.x, imag(uSol))

    savefig(plt, "directSolimag.pdf")


    println(a)
    println(b)


    # perform gmres-accelerated WaveHoltz
    #tol = 1e-14
    #u_0, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, tol)

    # perform WaveHoltz iterations:
    tol = 1e-9
    u_0, nIter = Waveholtz(simul, omega, fVals, alphas, g, tol)


    # compare with direct solution
    plt = heatmap(simul.y, simul.x, abs.(real(uSol) - u_0))
    savefig(plt, "realComparison.pdf")
    println(SEM_Wave_2d_Updated.LpNorm(simul, real(uSol), u_0, 2) / SEM_Wave_2d_Updated.LpNorm(simul, real(uSol), zeros(length(simul.y), length(simul.x)), 2))




    



end



function RnConvPlot()



    # domain
    R = 20.0
    xl = -R
    xr = R
    yl = -R
    yr = R


    N = 8   # degree of interpolation polynomials

    spaceTol = 1e-3



    omega = 5.1

    lambda = 1/(omega/(2*pi))                            # wavelength
    N_lambda = (xr-xl)/lambda                               # length of domain in wavelengths
    PPW = pi*(N_lambda / spaceTol)^(1/N)                    # points per wavelength
    dofs_x = ceil(PPW * N_lambda)                           # how many points do we want in space?

    K = Int(ceil((dofs_x - 1)/N))                           # the corresponding number of elements!


    maxIter = 250

    sVals = [1.5, 2, 3, 4, 5, 6]

    println("using $K by $K elements for frequency $omega")
    Kx = K # number of elements in x-direction
    Ky = K # number of elements in y-direction

    c_square(x, y) = 1

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    simul.exactStep = false

    nx = length(simul.x)
    ny = length(simul.y)

    # forcing
    S(x, y) = exp.(-((y .- (yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
    f(x, y) = omega^(5/2) * S(omega * x, omega * y)
    fVals = f.(simul.x, simul.y')

    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] #impedance all around

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    

    a = simul.stepCoeffs[1]
    b = simul.stepCoeffs[2]

    n = Int(length(simul.x) * length(simul.y))

    H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)

    uSolVec = H\F

    uDirSol = reshape(uSolVec, ny, nx)


    ############### WaveHoltz solve ###############
    #Waveholtz-specific parameters for the wave solver
    simul.fVals = -fVals
    simul.omega = omega
    simul.Tend = 2*pi/omega

    # WaveHoltz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))

    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))
    
    nsteps = Integer(2 * nsteps)
    
    # starting guess
    uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = zeros(length(simul.y), length(simul.x))

    # reset filtered solutions for WaveHoltz
    simul.uFiltered = zeros(length(simul.y), length(simul.x))
    simul.uDerFiltered = zeros(length(simul.y), length(simul.x))


    uSolDiffs = zeros(maxIter, length(sVals))

    for nIter = 1:maxIter

        SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g)
        
        uStart = simul.uFiltered
        #uStartDer = simul.uDerFiltered
        uStartDer = zeros(length(simul.y), length(simul.x))

        println("frequency: $omega || iteration: $nIter")

        for j = 1:length(sVals)

            s = sVals[j]
            bracket(x, y) = (1 + x^2 + y^2)^(-s/2)
            bracketMat = bracket.(simul.x, simul.y')
            
            uSolDiffs[nIter, j] = SEM_Wave_2d_Updated.LpNorm(simul, real(uDirSol)  .* bracketMat, uStart .* bracketMat, 2)

        end

    end




    s = 2.0
    bracket(x, y) = (1 + x^2 + y^2)^(-s/2)
    bracketMat = bracket.(simul.x, simul.y')


    plt = heatmap(simul.y, simul.x, bracketMat.*(real(uDirSol) - uStart))

    savefig(plt, "RealPartDiff.pdf")

    plt = heatmap(simul.y, simul.x, bracketMat.*(imag(uDirSol) - uStartDer))

    savefig(plt, "ImagPartDiff.pdf")


    plt = plot()

    for j = 1:length(sVals)

        plt = plot!(LinRange(1, maxIter, maxIter), uSolDiffs[:, j], xscale=:log10, yscale=:log10, label="s = $(sVals[j])")

    end

    plot!(LinRange(1, maxIter, maxIter), uSolDiffs[1, 1] ./ (sqrt.(LinRange(1, maxIter, maxIter))), xscale=:log10, yscale=:log10, label=L"O(n^{-1/2})", linestyle=:dash)



    savefig(plt, "ErrorsPlot.pdf")


    matwrite("RnConvData.mat", Dict(
    "H" => H,
    "F" => F,
    "sVals" => sVals,
    "U_whi" => simul.uNow,
    "uNorms" => uSolDiffs,
    "xVals" => simul.x,
    "yVals" => simul.y
    ))
    
end



function FrequencySweep()


    # domain
    R = 20.0
    xl = -R
    xr = R
    yl = -R
    yr = R

    N = 8   # degree of interpolation polynomials

    spaceTol = 1e-2


    maxIter = 200

    omegaMax = 12
    omegaCount = 4
    
    omegas = LinRange(3, omegaMax, omegaCount)

    sVals = [0.5, 1, 1.5, 2, 3, 4, 5, 6]                        # the different s parameters to use for norms

    waveholtzData = zeros(length(omegas), Integer(maxIter), length(sVals))

    dirSolNorms = zeros(length(omegas))


    for j = 1:omegaCount

        
        c_square(x, y) = 1.0                                           # constant wavespeed


        omega = omegas[j]
        lambda = 1/(omega/(2*pi))                               # wavelength
        N_lambda = (xr-xl)/lambda                               # length of domain in wavelengths

        N = 8                                                   # order of method

        PPW = pi*(N_lambda / spaceTol)^(1/N)                    # points per wavelength
        dofs_x = ceil(PPW * N_lambda)                           # how many points do we want in space?
        
        K = Int(ceil((dofs_x - 1)/N))                           # the corresponding number of elements!

        println("using $K by $K elements for frequency $omega")
        Kx = K                                                  # number of elements in x-direction
        Ky = K                                                  # number of elements in y-direction

        simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

        g = zeros(length(simul.y), length(simul.x))

        println("Using $K elements, i.e. $(length(simul.x)*length(simul.y)) points in space")

        S(x, y) = exp.(-((y .- (yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
        f(x, y) = omega^(5/2) * S(omega * x, omega * y)
        fVals = f.(simul.x, simul.y')                          # forcing

        alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)]  # impedance all around
        
        # construct the matrices:
        H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)

        uSolVec = H\F
        uDirSol = reshape(uSolVec, length(simul.y), length(simul.x))

        dirSolNorms[j] = SEM_Wave_2d_Updated.LpNorm(simul, real(uDirSol), zeros(length(simul.y), length(simul.x)), 2)




        ############### WaveHoltz solves ###############
        # Waveholtz-specific parameters for the wave solver
        simul.fVals = -fVals
        simul.omega = omega
        simul.Tend = 2*pi/omega
        

        # WaveHoltz-appropriate parameters for the wave solver
        Tend = 2*pi/omega
        delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
        delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
        cMax = sqrt(maximum(simul.c_square))

        nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))


        # starting guess
        uStart = zeros(length(simul.y), length(simul.x))
        uStartDer = zeros(length(simul.y), length(simul.x))

        # reset filtered solutions for WaveHoltz
        simul.uFiltered = zeros(length(simul.y), length(simul.x))
        simul.uDerFiltered = zeros(length(simul.y), length(simul.x))

        for nIter = 1:maxIter

            SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g)
            
            uStart = simul.uFiltered
            #uStartDer = simul.uDerFiltered
            uStartDer = zeros(length(simul.y), length(simul.x))
            println("frequency: $omega || iteration: $nIter")


            for k = 1:length(sVals)

                s = sVals[k]
                
                bracket(x, y) = (1 + x^2 + y^2)^(-s/2)

                bracketMat = bracket.(simul.x, simul.y')

                waveholtzData[j, nIter, k] = SEM_Wave_2d_Updated.LpNorm(simul, real(uDirSol)  .* bracketMat, uStart .* bracketMat, 2)

            end

        end

    end



    matwrite("NewFreqSweepData.mat", Dict(
    "omegas" => omegas,
    "sVals" => sVals,
    "uNorms" => dirSolNorms,
    "waveholtzData" => waveholtzData,
    ))





end


function HelmholtzMatrixTest()


    # domain
    R = 10.0
    xl = -R
    xr = R
    yl = -R
    yr = R


    N = 8   # degree of interpolation polynomials

    spaceTol = 1e-3



    omega = 5.1

    lambda = 1/(omega/(2*pi))                            # wavelength
    N_lambda = (xr-xl)/lambda                               # length of domain in wavelengths
    PPW = pi*(N_lambda / spaceTol)^(1/N)                    # points per wavelength
    dofs_x = ceil(PPW * N_lambda)                           # how many points do we want in space?

    K = Int(ceil((dofs_x - 1)/N))                           # the corresponding number of elements!


    maxIter = 250

    sVals = [1.5, 2, 3, 4, 5, 6]

    println("using $K by $K elements for frequency $omega")
    Kx = K # number of elements in x-direction
    Ky = K # number of elements in y-direction

    Kx = 20
    Ky = 20

    # uniform wave speed
    c_square(x, y) = 1

    # create SEM_Wave object to store data
    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    nx = length(simul.x)
    ny = length(simul.y)

    # forcing
    S(x, y) = exp.(-((y .- (yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
    f(x, y) = omega^(5/2) * S(omega * x, omega * y)
    fVals = f.(simul.x, simul.y')

    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] # impedance boundary conditions on all sides

    g = zeros(length(simul.y), length(simul.x)) # homogeneous impedance conditions
    simul.g = g
    
    to = TimerOutput()
    
    for j = 1:5

        # calculate the Helmholtz operators in two ways
        @timeit to "new" H1, F1 = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)
        @timeit to "old" H2, F2 = SEM_Wave_2d_Updated.HelmholtzMatrixOld(simul, omega, alphas, fVals)

    end

    show(to)

    #=
    plt = heatmap(abs.(H1-H2))
    savefig(plt, "test.png")

    println("Should be zero: $(norm(H1-H2, 2))") # is not zero!
    println("Should be zero: $(norm(F1-F2, 2))") # is zero!
    =#

end