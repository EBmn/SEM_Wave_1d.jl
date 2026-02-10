

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
    
    xl = -2.123/5
    xr = 4.0234/5
    yl = -3.234567/5
    yr = 1.5678/5

    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0




    N = 8   # degree of interpolation polynomials
    Kx = 8 # number of elements in x-direction
    Ky = 8 # number of elements in y-direction

    heaviside(x) = 0.5 * (sign(x) + 1)

    c_square(x, y) = 1 - 0.999*((heaviside(x-0.2) - heaviside(x-0.3)) .* (heaviside(y+0.9) - heaviside(y-0.9)) + (heaviside(x+0.3) - heaviside(x+0.2)) .* (heaviside(y+0.9) - heaviside(y-0.9)))



    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    Tend = 0.25
    #Tend = 0.034546
    nsteps = 200

    #=
    simul.useMMS = true
    simul.MMS_j.type = [3 3 3]

    simul.MMS_j.coeff[1, 1] = 1
    #simul.MMS_j.coeff[1, 2] = 3*pi * 2 / (xr - xl)
    simul.MMS_j.coeff[1, 2] = 0.0
    #simul.MMS_j.coeff[1, 3] = -3*pi*(xr + xl) / (xr - xl)
    #simul.MMS_j.coeff[1, 3] = 1.0
    simul.MMS_j.coeff[1, 3] = 3.0
    
    simul.MMS_j.coeff[2, 1] = 1
    #simul.MMS_j.coeff[2, 2] = 5*pi * 2 / (yr - yl)
    simul.MMS_j.coeff[2, 2] = 0.0
    #simul.MMS_j.coeff[2, 3] = -5*pi*(yr + yl) / (yr - yl)
    simul.MMS_j.coeff[2, 3] = 5.0
    #simul.MMS_j.coeff[2, 3] = 0.0

    #simul.MMS_j.coeff[2, 3] = 2.0
    #simul.MMS_j.coeff[2, 3] = 6.0

    simul.MMS_j.coeff[3, 1] = 1
    #simul.MMS_j.coeff[3, 2] = 1
    simul.MMS_j.coeff[3, 2] = 0
    #simul.MMS_j.coeff[3, 2] = 4*pi
    #simul.MMS_j.coeff[3, 3] = 0
    simul.MMS_j.coeff[3, 3] = 2


    simul.MMS_c.type = [3 3]

    simul.MMS_c.coeff[1, 1] = 1
    simul.MMS_c.coeff[1, 2] = 0.0
    simul.MMS_c.coeff[1, 3] = 1
    #simul.MMS_c.coeff[1, 1] = 1
    #simul.MMS_c.coeff[1, 2] = (xr+xl)/2
    #simul.MMS_c.coeff[1, 3] = 0.05


    simul.MMS_c.coeff[2, 1] = 1
    #simul.MMS_c.coeff[2, 2] = yl-1
    simul.MMS_c.coeff[2, 2] = 0
    simul.MMS_c.coeff[2, 3] = 1
    #simul.MMS_c.coeff[2, 1] = 1
    #simul.MMS_c.coeff[2, 2] = (yr+yl)/2
    #simul.MMS_c.coeff[2, 3] = 0.05


    =#


    #simul.useMMS = true
    simul.MMS_j.type = [3 3 3]

    simul.MMS_j.coeff[1, 1] = 1
    simul.MMS_j.coeff[1, 2] = 0.567
    simul.MMS_j.coeff[1, 3] = 2.0

    simul.MMS_j.coeff[2, 1] = 1
    simul.MMS_j.coeff[2, 2] = 0.0
    simul.MMS_j.coeff[2, 3] = 5.0

    simul.MMS_j.coeff[3, 1] = 1
    simul.MMS_j.coeff[3, 2] = 0
    simul.MMS_j.coeff[3, 3] = 1



    # troginometric test
    #=
    simul.MMS_j.type = [2 2 2]
    simul.MMS_j.coeff[1, 1] = 1.98765
    simul.MMS_j.coeff[1, 2] = 1 + 0.57721566490153286060651209008240243104215933593992
    simul.MMS_j.coeff[1, 3] = log(20)
    
    simul.MMS_j.coeff[2, 1] = 0.12345
    simul.MMS_j.coeff[2, 2] = sqrt(17)
    simul.MMS_j.coeff[2, 3] = 10*exp(-pi)

    simul.MMS_j.coeff[3, 1] = (1+sqrt(5))/2
    simul.MMS_j.coeff[3, 2] = 5*pi/4
    simul.MMS_j.coeff[3, 3] = 4
    =#


    
    
    #simul.MMS_j.coeff[2, 1] = 1
    #simul.MMS_j.coeff[2, 2] = 5*pi * 2 / (yr - yl)
    #simul.MMS_j.coeff[2, 2] = 0.02345
    #simul.MMS_j.coeff[2, 3] = -5*pi*(yr + yl) / (yr - yl)
    #simul.MMS_j.coeff[2, 3] = 1.0
    #simul.MMS_j.coeff[2, 3] = 6.0
    
    #simul.MMS_j.coeff[3, 1] = 1
    #simul.MMS_j.coeff[3, 2] = 1
    #simul.MMS_j.coeff[3, 3] = 3
    


    simul.MMS_c.type = [3 3]

    simul.MMS_c.coeff[1, 1] = 1
    #simul.MMS_c.coeff[1, 2] = 0
    simul.MMS_c.coeff[1, 2] = (xr+xl)/2
    simul.MMS_c.coeff[1, 3] = 2
    #simul.MMS_c.coeff[1, 1] = 1
    
    #simul.MMS_c.coeff[1, 3] = 0.05

    simul.MMS_c.coeff[2, 1] = 1
    simul.MMS_c.coeff[2, 2] = 0
    #simul.MMS_c.coeff[2, 1] = 1
    #simul.MMS_c.coeff[2, 2] = (yr+yl)/2
    simul.MMS_c.coeff[2, 3] = 2



    ###### set up the problem ######
    #fVals = cos.(3*pi*simul.x) + 5*exp.(-50*(simul.x .- (xr-xl)/2).^2) + 7*exp.(-70*(simul.x .- (xr-xl)/1.2).^2)
    #fVals = 1e2*exp.(-((simul.y .- 0.2*(yr + yl))/0.01).^2) * exp.(-((simul.x .- 0.7*(xr + xl))/0.01).^2)'
    fVals = 1e4*exp.(-((simul.y .- 0.2*(yr+yl))/0.03).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.03).^2)'
    #fVals = zeros(length(simul.y), length(simul.x))

    bc = [0.0, 1.0]
    a = 1/sqrt(2)
    #a = 0.0
    b = sqrt(1-a^2)
    bc = [a, b]

    g = zeros(length(simul.y), length(simul.x))


    simul.g = g
    alpha = simul.bc[1]
    beta = simul.bc[2]

    omega = 7.1

    Tend = 5.0
    nsteps = 2500

    # construct initial data
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    #uStart = exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'
    #uStart = 2 .- (abs.(simul.y) .+ abs.(simul.x'))

    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    #uStartDer = 1e3*exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'

    reference = zeros(length(simul.y), length(simul.x))

    if simul.useMMS
        for s = 1:length(simul.y)
            for t = 1:length(simul.x)
                uStart[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 0, simul.MMS_j)
                uStartDer[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 1, simul.MMS_j)
                reference[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], Tend, 0, 0, 0, simul.MMS_j)
            end
        end
    end

    useMMS = false
    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, true, 10, 5)

    #plt2 = surface(simul.x, simul.y[end:-1:1], log10.(abs.(simul.uNow-reference)))
    #savefig(plt2, "log10diff")
    #savefig(plt2, "log10diff")

    #logDiff = log10.(abs.(simul.uNow - reference))
    #plt = surface(simul.y, simul.x, logDiff)
    #plt = surface([1;2;3], [4;5;6], [1 0 2; 0 0 0; 3 0 4])
    #savefig(plt, "mmsLog10Diff")

    #gAppx, gMMS = SEM_Wave_2d.BoundaryMMS(simul, 1, bc)
    #mmsRef, appx = SEM_Wave_2d.InitialiseMMS(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)
    #plt = surface(simul.y, simul.x, log10.(abs.(gMMS-gAppx)))
    #plt = plot(simul.x, log10.(abs.(gMMS[1, :]-gAppx[1, :])))
    #plt = plot(simul.x, gMMS[1, :])
    
    #plt = plot(simul.x, gMMS[end, :])
    #plot!(simul.x, gAppx[end, :])
    #plt = plot(simul.y, gMMS[:, end])
    #plt = surface(simul.x, simul.y[end:-1:1], log10.(abs.(mmsRef - appx)))
    #savefig(plt, "boundaryMMS_tester")

    #plt2 = surface(simul.x, simul.y[end:-1:1], uStart)
    #savefig(plt2, "tZeroRef")

    #println(maximum(gAppx))
    #println(minimum(gAppx))

    #println(maximum(gMMS))

    #println(minimum(gMMS))

    #appx, mmsRef = SEM_Wave_2d.LaplaceMMS(simul)
    #logDiff = log10.(abs.(mmsRef-appx))

    #println((mmsRef-appx)[2:end-1, 2:end-1])
    #println(mmsRef)
    #println(appx)

    #plt = surface(simul.x[2:end-1], simul.y[end-1:-1:2], logDiff[2:end-1, 2:end-1])
    #plt = surface(simul.x[2:end-1], simul.y[end-1:-1:2], (mmsRef-appx)[2:end-1, 2:end-1])
    #plt = surface(simul.x[2:end-1], simul.y[2:end-1], appx[2:end-1, 2:end-1])
    #plt = surface(simul.x[2:end-1], simul.y[2:end-1], mmsRef[2:end-1, 2:end-1])
    #plt = surface(simul.x, simul.y[end:-1:1], simul.c_square)
    #println(maximum(abs.(appx)))
    #println(maximum(abs.(mmsRef)))
    #plt = surface(simul.x, simul.y[end:-1:1], log10.(abs.(appx)))
    
    #println(logDiff[1, 1])
    #println(logDiff[1, end])
    #println(logDiff[end, 1])
    #println(logDiff[end, end])

    #savefig(plt, "Lap_tester")

    #println(simul.QuadPoints)

    #spy(simul.M_b)

    # tests the boundary integral evaluated using simul.M_b
    #=
    p = 15 
    q = 5
    vals = (simul.x'.^p .+ simul.y[end:-1:1].^q)     # integrand
    appxInt = sum(vals .* simul.M_b)                 # approximate line integral
    exactInt = (xr^(p+1)-xl^(p+1))*2/(p+1) + (xl^p + xr^p)*(yr - yl) + (yr^(q+1)-yl^(q+1))*2/(q+1) + (yl^q + yr^q)*(xr - xl) # exact line integral
    println("Approximate line integral, " * string(appxInt) * " should be close to the exact value, " * string(exactInt))
    =#
    
    #println((xr^(p+1)-xl^(p+1))*2/(p+1) + (xl^p + xr^p)*(yr - yl) + (yr^(q+1)-yl^(q+1))*2/(q+1) + (yl^q + yr^q)*(xr - xl))  
    
    
    #plt = surface(simul.y, simul.x, simul.M_b)
    #savefig(plt, "M_b_tester")
    #println(simul.QuadWeights)

    

end

