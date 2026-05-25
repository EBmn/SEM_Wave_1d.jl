



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


    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0



    N = 8   # degree of interpolation polynomials
    Kx = 60 # number of elements in x-direction
    Ky = 60 # number of elements in y-direction

    c_square(x, y) = 1.0

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    ###### set up the problem ######
    ny = 20
    nx = 20
    yExpRange = LinRange(0, ny, ny+1)
    xExpRange = LinRange(0, nx, nx+1)

    poinConstLow = -Inf
    normQuotients = zeros(ny+1, nx+1)

    for i = 1:length(xExpRange)
        for j = 1:length(yExpRange)

            xExp = Int(round(2*xExpRange[i]))
            yExp = Int(round(2*yExpRange[j] + 1))

            u = (simul.y).^(yExp) * (simul.x').^xExp # this function has zero mean over our domain, so \|u\|_H^1 \leq C \|\nabla\|_L^2, where the latter can be calculated with GradIntegral

            u_L2 = SEM_Wave_2d.LpNorm(simul, u, zeros(length(simul.y), length(simul.x)), 2)
            u_H1 = SEM_Wave_2d.SobolevNorm(simul, u)
            uGradNorm = (SEM_Wave_2d.GradIntegral(simul, u, u))^(1/2)

            quotient = u_H1 / uGradNorm
            normQuotients[j, i] = quotient

            if (quotient > poinConstLow)
                poinConstLow = quotient
            end


            analyticGradL2 = sqrt((16*xExpRange[i]^2) / ((4*xExpRange[i] - 1)*(4*yExpRange[j] + 3)) + (4*(2*yExpRange[j] + 1)^2 )/ ((4*xExpRange[i] + 1)*(4*yExpRange[j] + 1)))
            analyticL2 = sqrt(4 / ((4*xExpRange[i] + 1)*(4*yExpRange[j] + 3)))
            analyticH1 = analyticL2 + analyticGradL2

            analyticQuotient = analyticH1/analyticGradL2

            println("y^" * string(yExp) * "x^" * string(xExp) * " gives quotient " * string(quotient))
            println("corresponding analytical quotient: " * string(analyticQuotient) * ", reldiff of " * string(abs((quotient-analyticQuotient)/analyticQuotient)))

            #=
            println("H1 relative error: " * string(abs((u_H1 - analyticH1)/analyticH1)))
            println("u L2 relative error: " * string(abs((u_L2 - analyticL2)/analyticL2)))
            println("uGrad relative L2 error: " * string(abs((uGradNorm - analyticGradL2))/analyticGradL2))
            =#



        end
    end

    println("empirical lower bound of Poincare constant: " *string(poinConstLow))


    plt = heatmap(2*yExpRange .+ 1, 2*xExpRange, normQuotients)

    savefig(plt, "normQuotientHeatMap")

end

