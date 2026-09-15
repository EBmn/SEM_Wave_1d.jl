


using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using TimerOutputs
using Base.Threads
using MAT

include("SEM_Wave_2d_Updated.jl")
using .SEM_Wave_2d_Updated



include("MMS.jl")
using .MMS


function normsTester()

    ###### setup ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    alphas = [0.0; 0.0; 0.0; 0.0]
    
    omega = 5.0

    N = 6   # degree of interpolation polynomials

    c_square_1(x, y) = 1 # first uniform wavespeed: GradIntegral and GradIntegralSeminorm should agree in this case


    Kx = 4
    Ky = 4

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square_1)
    wd = 0.2
    f = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)' # some function values. 
    g = 3*sin.(simul.y .- 0.1) * cos.(3*simul.x .+ 0.4)'


    gradIntegral1 = SEM_Wave_2d_Updated.GradIntegral(simul, f, g)
    gradIntegral2 = SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, f, g)

    sobNorm1 = SEM_Wave_2d_Updated.SobolevNorm(simul, g, 1.0)
    semNorm1 = SEM_Wave_2d_Updated.SeminormWH(simul, g, g)

    
    println("should be zero: $(abs(gradIntegral1 - gradIntegral2))")
    println("should be zero: $(abs(sobNorm1 - semNorm1))")
    

    c_square(x, y) = 1 + 0.7*sin.(x - sqrt(2))*sin.(y) # variable wavespeed: GradIntegral and GradIntegralSeminorm should no longer agree
    
    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    gradIntegral3 = SEM_Wave_2d_Updated.GradIntegral(simul, f, g)
    gradIntegral4 = SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, f, g)

    sobNorm2 = SEM_Wave_2d_Updated.SobolevNorm(simul, g, 1.0)
    semNorm2 = SEM_Wave_2d_Updated.SeminormWH(simul, g, g)

    println("should be zero: $(abs(gradIntegral1 - gradIntegral3))") # this function should not see the c^2, so that shouldn't matter
    println("should not be zero: $(abs(gradIntegral2 - gradIntegral4))") # this function should see the c^2, so these should not be equal.

    println("should be zero: $(abs(sobNorm1 - sobNorm2))")
    println("should not be zero: $(abs(semNorm1 - semNorm2))")


    f = (1.0 .+ simul.x'.^2 .+ simul.y.^2)
    g = (simul.x'.^2 .- simul.y.^2)
    gradInt = SEM_Wave_2d_Updated.GradIntegral(simul, f, g)
    println("should be zero: $(gradInt)")
    gradInte = SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, f, g)
    println("should probably not be zero: $(gradInte)")


    c_square_2(x, y) = 4 + x + y
    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square_2)
    gradInteg = SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, f, g)
    println("should be zero: $(gradInteg)")


end
