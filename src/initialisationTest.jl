


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
    omega = 9.225806451612902

    c_square(x, y) = 1
    c_min = 1.0 ### has to be set properly! ###

    L = maximum([xr-xl, yr-yl]) # length of the domain
    
    dofs_x = 80 # how many points do we want per dimension?
    numberOfElements = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!


    Kx = numberOfElements # number of elements in x-direction
    Ky = numberOfElements # number of elements in y-direction

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    fVals = ones(length(simul.y), 1) * cos.(3*pi*simul.x)'
    
    g = zeros(length(simul.y), length(simul.x))        
 

    ### from SEM_Wave_2d.Waveholtz:
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    spaceStep = minimum([delta_x delta_y])
    nsteps = Integer(ceil(1.5 * Tend * (1/delta_x + 1/delta_y)))

    uStart = ones(length(simul.y), 1) * cos.(4*pi*simul.x)'
    #uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = ones(length(simul.y), 1) * cos.(5*pi*simul.x)'
    #uStartDer = zeros(length(simul.y), length(simul.x))

    SEM_Wave_2d.Initialise!(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)

    #= 
    this means that the solution to the wave equation with zero initial condition is precisely 
    (1/(omega^2 - 9*pi^2)) * (cos(omega*t) - cos(3*pi*t)) * cos(3*pi*x) + cos(4*pi*t)*cos(4*pi*x) + (1/(5*pi)) * sin(5*pi*t)*cos(5*pi*x)
    =#

    t = -Tend/nsteps
    println(t)
    uPrevTrue = (1/(omega^2 - 9*pi^2)) * (cos(omega*t) - cos(3*pi*t)) * ones(length(simul.y), 1) * cos.(3*pi*simul.x)' + cos(4*pi*t) * ones(length(simul.y), 1) * cos.(4*pi*simul.x)' + 
                (1/(5*pi)) * sin(5*pi*t) * ones(length(simul.y), 1) * cos.(5*pi*simul.x)'

    println(maximum(simul.uPrev-uPrevTrue))

    plt = surface(simul.x, simul.y[end:-1:1], simul.uPrev-uPrevTrue)
    #plt = surface(simul.x, simul.y[end:-1:1], uPrevTrue)
    #plt = surface(simul.x, simul.y[end:-1:1], simul.uPrev)
    savefig(plt, "initTest.png")

end