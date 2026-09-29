
using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using Base.Threads
using Dates
using SparseArrays
using MAT


include("SEM_Wave_2d_Updated.jl")
using .SEM_Wave_2d_Updated


include("MMS.jl")
using .MMS




function PlotData()

    ###### setup ######      
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 2.0


    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 0.0]
    
    omega = 75.0
    tol = 1e-9

    N = 10   # degree of interpolation polynomials


    x0 = 0.0; y0 = 0.75; x1 = 0.25; y1 = 0.75; x2 = 0.5; y2 = 0.75; x3 = 0.75; y3 = 0.75; x4 = 1.0; y4 = 0.75;
    c_square(x, y) = 1 - 0.75*Float64(((x-x0)/0.9)^2 + ((y-y0)/0.5)^2 < 0.1^2) - 0.75*Float64(((x-x1)/0.9)^2 + ((y-y1)/0.5)^2 < 0.1^2) - 0.75*Float64(((x-x2)/0.9)^2 + ((y-y2)/0.5)^2 < 0.1^2) - 0.75*Float64(((x-x3)/0.9)^2 + ((y-y3)/0.5)^2 < 0.1^2)  - 0.75*Float64(((x-x4)/0.9)^2 + ((y-y4)/0.5)^2 < 0.1^2)

    Kx = 20
    Ky = 40


    # create simul object and forcing
    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    
    g = zeros(length(simul.y), length(simul.x))       
    w = 0.05
    fVals = (omega^2) * exp.(-((simul.y .- (yr-yl)/10)./w).^2) * exp.(-((simul.x .- (xr-xl)/2)./w).^2)'

    u_re, u_im, hist = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, tol)
    u_im = u_im/omega

    # store some stuff for later
    c_data = simul.c_square
    xVals = simul.x
    yVals = simul.y

    
    matwrite("ImpedanceForFun.mat", Dict(
    "u_re" => u_re,
    "u_im" => u_im,
    "x" => xVals,
    "y" => yVals,
    "f" => fVals,
    "c_data" => c_data,
    "omega" => omega,
    "N" => N
    ))

end