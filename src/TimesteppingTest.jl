
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

    omega = 7.1
    #omega = 2.1


    N = 8   # degree of interpolation polynomials
    Kx = 12 # number of elements in x-direction
    Ky = 12 # number of elements in y-direction

    heaviside(x) = 0.5 * (sign(x) + 1)
    box(x, a, b) = heaviside(x-a) - heaviside(x-b)

    c_square(x, y) = 1 - 0.9*box(x, 0.2, 0.5)*box(y, 0.3, 0.6) - 0.9*box(x, -0.5, -0.2)*box(y, -0.6, -0.3)

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    nx = length(simul.x)
    ny = length(simul.y)
    
    f(x, y) = exp.(-((y .- 0.2*(yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
    #f(x, y) = omega^2

    fVals = f.(simul.x, simul.y')

    #alphas = [0.0; 0.0; 0.0; 0.0] # Neumann all around
    #alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] # impedance all around
    alphas = [1/sqrt(2); 0.0; 1/sqrt(2); 0.1] # mixed conditions
    

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


#=
    # perform gmres-accelerated WaveHoltz iterations
    tol = 1e-14

    simul.exactStep = false
    v_0, v_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, tol)

    simul.exactStep = true
    u_0, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, tol)
=#

    



       
    # perform WaveHoltz iterations:
    tol = 5e-15

    simul.exactStep = false
    v_0, v_1, nIter1 = SEM_Wave_2d_Updated.Waveholtz(simul, omega, fVals, alphas, g, tol)

    simul.exactStep = true
    u_0, u_1, nIter2 = SEM_Wave_2d_Updated.Waveholtz(simul, omega, fVals, alphas, g, tol)



    # compare with direct solution
    plt = heatmap(simul.y, simul.x, log10.(abs.(real(uSol) - u_0)))
    savefig(plt, "realComparison1.pdf")
    #println("real(uSol) - u_0: $(norm(real(uSol) - u_0))")
    println("should be zero: $(norm(real(uSol) - u_0))")
    
    #println(norm(real(uSol) - u_0) / norm(u_0))
    
    plt = heatmap(simul.y, simul.x, log10.(abs.(real(uSol) - v_0)))
    savefig(plt, "realComparison2.pdf")
    #println("real(uSol) - v_0: $(norm(real(uSol) - v_0))")
    println("should not be zero: $(norm(real(uSol) - v_0))")
    #println(SEM_Wave_2d_Updated.LpNorm(simul, real(uSol), u_0, 2) / SEM_Wave_2d_Updated.LpNorm(simul, real(uSol), zeros(length(simul.y), length(simul.x)), 2))

    #plt = heatmap(simul.y, simul.x, log10.(abs.(v_0 - u_0)))
    #savefig(plt, "realComparison3.pdf")
    #println("v_0 - u_0: $(norm(v_0 - u_0))")
    
    plt = heatmap(simul.y, simul.x, log10.(abs.(imag(uSol) - u_1/omega)))
    savefig(plt, "imagComparison1.pdf")
    #println("imag(uSol) - u_1/omega: $(norm(imag(uSol) - u_1/omega))")
    println("should be zero: $(norm(imag(uSol) - u_1/omega))")
    
    plt = heatmap(simul.y, simul.x, log10.(abs.(imag(uSol) - v_1/omega)))
    savefig(plt, "imagComparison2.pdf")
    #println("imag(uSol) - v_1/omega: $(norm(imag(uSol) - v_1/omega))")
    println("should not be zero: $(norm(imag(uSol) - v_1/omega))")
    
    #plt = heatmap(simul.y, simul.x, log10.(abs.(v_1 - u_1)))
    #savefig(plt, "imagComparison3.pdf")
    #println("v_1 - u_1: $(norm(v_1 - u_1))")

    
    #println(norm(imag(uSol) - u_1/omega) / norm(imag(uSol)))
    

    #println(norm(uSol - (u_0 + 1im*u_1/omega)))

    plt = heatmap(simul.y, simul.x, real(uSol))
    savefig(plt, "tester1.pdf")

    plt = heatmap(simul.y, simul.x, u_0)
    savefig(plt, "tester2.pdf")

    plt = heatmap(simul.y, simul.x, imag(uSol))
    savefig(plt, "tester3.pdf")

    plt = heatmap(simul.y, simul.x, u_1/omega)
    savefig(plt, "tester4.pdf")

end


function test1()

    # tests that there is a difference at all between the old timestepping scheme (stepcoeffs = 1) and the new one.
    # does this by running the same WaveholtzGMRES call twice, but with different timestepping schemes, and then comparing the results
    
    # domain
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0

    omega = 5.1


    N = 3   # degree of interpolation polynomials
    Kx = 4 # number of elements in x-direction
    Ky = 4 # number of elements in y-direction

    c_square(x, y) = 1

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    nx = length(simul.x)
    ny = length(simul.y)
    
    f(x, y) = exp.(-((y .- 0.2*(yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
    fVals = f.(simul.x, simul.y')

    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] #impedance all around

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    

    



    # perform WaveHoltz iterations: modified timestepping
    println(simul.stepCoeffs[1])
    println(simul.stepCoeffs[2])


    tol = 1e-14
    u_0, u_1, history0 = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, tol)


    # perform WaveHoltz iterations: original leapfrog

    simul.exactStep = false

    println(simul.stepCoeffs[1])
    println(simul.stepCoeffs[2])
    v_0, v_1, history1 = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, tol)
    println(simul.stepCoeffs[1])
    println(simul.stepCoeffs[2])


    println("should be nonzero: " * string(norm(u_0 - v_0)))
    println("should be nonzero: " * string(norm(u_1 - v_1)))

end


function test2()

    # tests the exact timestepping on a manufactured function for which the timestepping should be exact
    # neumann conditions, so for real forcings the solution should be real
    # we find the Helmholtz solution V using a direct solver; timestepping from w(0) = V, w_t(0) = 0 should produce the solution Vcos(omega*t) to machine precision

    # domain
    scale = 0.05
    xl = -scale
    xr = scale
    yl = -scale
    yr = scale

    omega = 10/scale


    N = 8   # degree of interpolation polynomials
    Kx = 12 # number of elements in x-direction
    Ky = 12 # number of elements in y-direction

    c_square(x, y) = 1 + (x^2 + y^2)/scale^2

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    nx = length(simul.x)
    ny = length(simul.y)
    
    f(x, y) = exp.(-((y .- 0.2*(yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
    fVals = f.(simul.x, simul.y')

    alphas = [0.0; 0.0; 0.0; 0.0] # Neumann all around
    #alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] # impedance all around

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    
    ##### find Helmholtz solution directly #####
    H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)
    V = reshape(H\F, ny, nx) # the solution


    # test instead with an eigenfunction
    #V = cos.(pi*simul.x)' .* ones(length(simul.y), 1)
    #omega = 1.0*pi

    ##### perform wave solve with V as initial data #####
    uStart = real(V)
    uStartDer = omega*imag(V)

    Tend = 2*scale
    nsteps = 2400

    #simul.exactStep = false

    SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)


    ##### compare result with what we expect to see #####
    u_expected = real(V*exp(-im*omega*Tend))
    println("should be zero: " * string(norm(simul.uNow - u_expected)))
    
    plt = heatmap(simul.y, simul.x, abs.(simul.uNow - u_expected))
    savefig(plt, "timesteppingDiff.pdf")

end


function test3()

    # tests whether the filtering is exact

    # domain
    scale = 1.0
    xl = -scale
    xr = scale
    yl = -scale
    yr = scale

    omega = 7.1


    N = 8   # degree of interpolation polynomials
    Kx = 12 # number of elements in x-direction
    Ky = 12 # number of elements in y-direction

    c_square(x, y) = 1 + (x^2 + y^2)/scale^2

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    nx = length(simul.x)
    ny = length(simul.y)
    
    f(x, y) = exp.(-((y .- 0.2*(yr+yl))/0.1).^2) .* exp.(-((x .- (xr+xl))/0.1).^2)
    #f(x, y) = omega^2
    fVals = f.(simul.x, simul.y')

    #alphas = [0.0; 0.0; 0.0; 0.0] # Neumann all around
    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] # impedance all around

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    
    ##### find Helmholtz solution directly #####
    H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)
    V = reshape(H\F, ny, nx) # the solution


    # test instead with an eigenfunction
    #V = cos.(pi*simul.x)' .* ones(length(simul.y), 1)
    #omega = 1.0*pi

    ##### perform wave solve with V as initial data #####
    uStart = real(V)
    uStartDer = omega*imag(V)

    Tend = 2*pi/omega
    nsteps = 600

    # take one more step so that we filter the entirety of [0, 2*pi/omega]
    #dt = Tend/nsteps
    #Tend = Tend + dt
    #nsteps = nsteps + 1

    #simul.exactStep = false

    # with these initial conditions we expect simul.uFiltered and simul.uDerFiltered to equal real(V) and omega*imag(V) respectively
    SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)

    

    ##### compare result with what we expect to see #####
    u_expected = real(V)
    uDer_expected = omega*imag(V)
    println("maximum abs of imag(V): " * string(maximum(abs.(imag(V)))))
    println("maximum of V: " * string(maximum(u_expected)))
    println("difference between max(V) and min(V): " * string(maximum(u_expected) - minimum(u_expected)))



    println("should be zero (real part): " * string(norm(simul.uFiltered - u_expected)))
    println("should be zero pointwise (real part): " * string(simul.uFiltered[14, 15] - u_expected[14, 15]))
    println("should be zero (imag part): " * string(norm(simul.uDerFiltered - uDer_expected)))
    
    plt = heatmap(simul.y, simul.x, abs.(simul.uFiltered - u_expected))
    plt = plot(simul.y, simul.x, log10.(abs.(simul.uFiltered - u_expected)))
    savefig(plt, "filterError1.pdf")
    
    plt = heatmap(simul.y, simul.x, abs.(simul.uDerFiltered - uDer_expected))
    plt = plot(simul.y, simul.x, log10.(abs.(simul.uDerFiltered - uDer_expected)))
    savefig(plt, "filterError2.pdf")

end