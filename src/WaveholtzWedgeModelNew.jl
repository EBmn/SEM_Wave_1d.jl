


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




function frequencySweep()

    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0


    alphas = [0.0; 0.0; 0.0; 0.0]

    
    spaceTol = 1e-2
    systemTol = 1e-3

    #systemTol = 1e3


    #set up range of omegas
    
    #= 
    # setup that crashed
    spaceTol = 1e-6
    systemTol = 1e-3
    omegaSmall = 15.3225
    omegaBig = 20.0
    numberOfOmegas = 30
    =#

    
    h = 0.05
    omegaSmall = 21.8
    omegaBig = 22.0
    
    numberOfOmegas = Int(round((omegaBig-omegaSmall) / h + 1))
    omegas = collect(LinRange(omegaSmall, omegaBig, numberOfOmegas))


    #run for these values of omegas, record the number of iterations
    gmresIters = zeros(length(omegas))
    whIters = zeros(length(omegas))


    # GPT-generated progress tracking
    #nT = Threads.nthreads()
    nT = Threads.maxthreadid()
    current_omega = fill(NaN, nT)

    n = length(omegas)

    ch = Channel{Int}(n)

    # fill channel with work
    for j in eachindex(omegas)
        put!(ch, j)
    end

    close(ch)



    completed = falses(n)
    done = Ref(false)

    # start monitor
    task_monitor = Threads.@spawn monitor_progress(completed, n, done, current_omega)

    #tasks = Vector{Task}(undef, length(omegas))

    tasks = Vector{Task}(undef, Threads.nthreads()-1) # leave one thread for the monitor!

    KxVals = zeros(length(omegas))

    for t in 1:length(tasks)

        #println("about to spawn task $t")

        tasks[t] = Threads.@spawn begin

            
            try
                for j in ch

                    
                    tid = threadid()
                    omega = omegas[j]
                    
                    #println("task $t considering value $omega")

                    current_omega[tid] = omega

                    N = 6   # degree of interpolation polynomials

                    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.
                    c_min = 1000 # has to be chosen appropriately

                    # there are three parts of the domain
                    lambda_1 = 2100/(2*pi*omega)        # region 1: wavespeed of 2100 and max length of 500 units in y-direction
                    lambda_1 = 2100/(omega/(2*pi))
                    N1_lambda = 500 / lambda_1

                    lambda_2 = 1000/(2*pi*omega)        # region 2: wavespeed of 1000 and max length of 400 units in y-direction
                    lambda_2 = 1000/(omega/(2*pi))
                    N2_lambda = 400 / lambda_2

                    lambda_3 = 2900 / (2*pi*omega)      # region 3: wavespeed of 2900 and max length of 400 units in y-direction
                    lambda_3 = 2900 / (omega/(2*pi))
                    N3_lambda = 400 / lambda_3

                    N_lambda = N1_lambda + N2_lambda + N3_lambda # upper bound on domain length in wavelengths


                    PPW = pi*(N_lambda/spaceTol)^(1/N)
                    
                    dofs_x = ceil(PPW * N_lambda)       # how many points do we want per dimension? points = point/wavelength * wavelength

                    Kx = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!
                    Ky = Integer(ceil(Kx * (5/3))) # number of elements in y-direction


                    if Kx < 3
                        Kx = 3
                        Ky = Integer(ceil(Kx * (5/3)))
                    end

                    #=
                    if Kx < 5
                        Kx = 5
                        Ky = Integer(ceil(Kx * (5/3)))
                    end
                    =#
                    #println("$omega: $Kx, $Ky")

                    # There seems to be a difference in iteration counts for odd and even element counts. 
                    # This makes sure we stay in the same category (odd number of elements in each direction) throughout
                    if iseven(Kx)
                        Kx = Kx + 1
                    end
                    
                    if iseven(Ky)
                        Ky = Ky + 1
                    end


                    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
                    
                    fVals = 100 * exp.(-5*((simul.y.-123)./1000).^2) * exp.(-5*((simul.x .- 245.0)./600).^2)'

                    
                    
                    g = zeros(length(simul.y), length(simul.x))        

                    
                    

                    println("$omega: $Kx by $Ky")

                    
                    sol1, sol2, nIter = SEM_Wave_2d_Updated.Waveholtz(simul, omega, fVals, alphas, g, systemTol)
                    u_0, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, systemTol)


                    gmresIter = length(history.data[:resnorm])
                    

                    gmresIters[j] = gmresIter
                    whIters[j] = nIter
                    KxVals[j] = Kx

                    completed[j] = true
                    current_omega[tid] = NaN
                
                end
                    
            
            catch e
                println("Worker thread $(threadid()) crashed!")
                showerror(stdout, e)
                println()
                rethrow(e)
            end
            
        end

    end


    #fetch.(tasks)
    wait.(tasks)

    done[] = true
    
    wait(task_monitor)
    
    
    # plot the result

    plt = plot(omegas, whIters, yscale=:log10, label="Waveholtz iteration", legend=:bottomleft)
    plot!(omegas, gmresIters, yscale=:log10, label="gmres iterations", legend=:bottomleft)

    savefig(plt, "WedgeIterPlot")

    println("values of omega: ")
    println(omegas)
    println("iterations for convergence: whi")
    println(whIters)
    println("iterations for convergence: gmres")
    println(gmresIters)
    println("elements in x-direction:")
    println(KxVals)

end



function monitor_progress(completed, total, done, current_omega) # GPT-generated

    t_start = time()

    while !done[]
        sleep(10.0)

        n_done = count(completed)
        progress = n_done / total

        elapsed = time() - t_start
        rate = n_done / max(elapsed, 1e-6)   # jobs per second

        remaining = total - n_done
        eta = remaining / max(rate, 1e-6)

        # progress bar
        bar_length = 30
        filled = round(Int, bar_length * progress)
        bar = "█"^filled * "-"^(bar_length - filled)

        # clear screen
        print("\033[H\033[J")



        println("=== Simulation Progress ===")
        println("[$bar] $(round(progress*100, digits=1)) %")
        println("$n_done / $total completed")

        println("\nSpeed: $(round(rate, digits=2)) runs/sec")
        println("Elapsed: $(round(elapsed, digits=1)) s")
        println("ETA: $(round(eta, digits=1)) s")



        println("\nActive threads:")
        
        for tid in eachindex(current_omega)
            omega = current_omega[tid]
            if isnan(omega)
                #println("Thread $tid: idle")
            else
                println("Thread $tid: omega = $omega")
            end
        end


    end
end



function convHistory()

    ###### setup ######    
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0

    alphas = [0.0; 0.0; 0.0; 0.0]
    
    omega = 22.0

    spaceTol = 1e-3
    systemTol = 1e-12

    nIter = 200

    N = 6   # degree of interpolation polynomials

    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.

    Kx = 5
    Ky = 7



    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
    fVals = 100 * exp.(-5*((simul.y.-123)./1000).^2) * exp.(-5*((simul.x .- 245.0)./600).^2)'
    g = zeros(length(simul.y), length(simul.x))       

    ref_0, ref_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, systemTol)
    ref_1 = ref_1 / omega # so as to make ref_1 an approximation of Im(u)

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))
    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))

    uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = zeros(length(simul.y), length(simul.x))

    diffs_0 = zeros(nIter, 1)
    diffs_1 = zeros(nIter, 1)
    diffs = zeros(nIter, 1)

    for j = 1:nIter

        # perform WHI and measure H_1-differences with reference
        SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

        diffs_0[j] = SEM_Wave_2d_Updated.SobolevNorm(simul, simul.uFiltered - ref_0, omega)
        diffs_1[j] = SEM_Wave_2d_Updated.SobolevNorm(simul, simul.uDerFiltered - omega * ref_1, omega)
        diffs[j] = SEM_Wave_2d_Updated.SeminormWH(simul, simul.uFiltered - ref_0, simul.uDerFiltered - ref_1)


        uStart .= simul.uFiltered
        uStartDer .= simul.uDerFiltered 

        println("iter $j of $nIter")
        
    end
                        
    # plot the result

    plt = plot(LinRange(1, nIter, nIter), diffs_0.^2 + diffs_1.^2, yscale=:log10, label="total iteration error", legend=:bottomleft)
    plot!(LinRange(1, nIter, nIter), diffs_0, yscale=:log10, label="v_0 error")
    plot!(LinRange(1, nIter, nIter), diffs_1, yscale=:log10, label="v_1 error")
    plot!(LinRange(1, nIter, nIter), diffs, yscale=:log10, label="seminorm error")



    savefig(plt, "WedgeConvHistory")

    plt = heatmap(simul.x, simul.y, ref_0)
    savefig(plt, "wedgeRealPart.pdf")

    println("v_0 error: ")
    println(diffs_0)
    println("v_1 error: ")
    println(diffs_1)
    println("seminorm error: ")
    println(diffs)


end



function convHistoryWedge()

    ###### setup ######    
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0

    alphas = [0.0; 0.0; 0.0; 0.0]
    
    omega = 20.0

    nIter = 500

    N = 6   # degree of interpolation polynomials

    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.

    KxBase = 24
    KyBase = 40

    
    # the data we want: second column has the data from the finer spatial discretisation
    e_0_L2 = zeros(nIter, 2)
    e_1_L2 = zeros(nIter, 2)
    grad_e_0_L2 = zeros(nIter, 2)
    grad_e_1_L2 = zeros(nIter, 2)
    
    grad_c_e_0_L2 = zeros(nIter, 2)
    grad_c_e_1_L2 = zeros(nIter, 2)

    # the same thing but with inexact timestepping
    e_0_L2_inexact = zeros(nIter, 2)
    e_1_L2_inexact = zeros(nIter, 2)
    grad_e_0_L2_inexact = zeros(nIter, 2)
    grad_e_1_L2_inexact = zeros(nIter, 2)

    grad_c_e_0_L2_inexact = zeros(nIter, 2)
    grad_c_e_1_L2_inexact = zeros(nIter, 2)



    number_of_refinements = 2

    first_nx = (N-1)*KxBase + 1
    first_ny = (N-1)*KyBase + 1 

    u_re_coarse = zeros(first_nx, first_ny)
    u_im_coarse = zeros(first_nx, first_ny)
    xVals_coarse = zeros(first_nx, 1)
    yVals_coarse = zeros(first_ny, 1)
    fVals_coarse = zeros(first_nx, first_ny)

    second_nx = (N-1)*2*KxBase + 1
    second_ny = (N-1)*2*KyBase + 1

    u_re_fine = zeros(second_nx, second_ny)
    u_im_fine = zeros(second_nx, second_ny)
    xVals_fine = zeros(second_nx, 1)
    yVals_fine = zeros(second_ny, 1)
    fVals_fine = zeros(second_nx, second_ny)
    c_data_fine = zeros(second_nx, second_ny)
    
    u_re_H1_coarse = 0.0
    u_re_H1_fine = 0.0
    u_im_H1_coarse = 0.0
    u_im_H1_fine = 0.0
    u_seminorm_coarse = 0.0
    u_seminorm_fine = 0.0

    for j = 1:number_of_refinements


        # choose number of elements
        Kx = Integer(KxBase*j)
        Ky = Integer(KyBase*j)


        # create simul object and forcing
        simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
        
        #wd = 0.2
        #fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'
        fVals = 100 * exp.(-5*((simul.y.-123)./1000).^2) * exp.(-5*((simul.x .- 245.0)./600).^2)'

        g = zeros(length(simul.y), length(simul.x))       

        # reference solution
        H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)
        uSolVec = H\F
        uSol = reshape(uSolVec, length(simul.y), length(simul.x))
        u_re = real(uSol)
        u_im = imag(uSol)

        # store the reference solutions for later
        if (j==1)
            u_re_coarse = u_re
            u_im_coarse = u_im

            xVals_coarse = simul.x
            yVals_coarse = simul.y
            fVals_coarse = fVals

            u_re_H1_coarse = SEM_Wave_2d_Updated.SobolevNorm(simul, u_re, 1.0)
            u_im_H1_coarse = SEM_Wave_2d_Updated.SobolevNorm(simul, u_im, 1.0)
            u_seminorm_coarse = SEM_Wave_2d_Updated.SeminormWH(simul, u_re, u_im)
    
        else
            u_re_fine = u_re
            u_im_fine = u_im

            xVals_fine = simul.x
            yVals_fine = simul.y
            fVals_fine = fVals

            u_re_H1_fine = SEM_Wave_2d_Updated.SobolevNorm(simul, u_re, 1.0)
            u_im_H1_fine = SEM_Wave_2d_Updated.SobolevNorm(simul, u_im, 1.0)
            u_seminorm_fine = SEM_Wave_2d_Updated.SeminormWH(simul, u_re, u_im)

            c_data_fine = simul.c_square
            
        end

        


        # Waveholtz-appropriate parameters for the wave solver
        Tend = 2*pi/omega
        delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
        delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
        cMax = sqrt(maximum(simul.c_square))
        nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))

        uStart = zeros(length(simul.y), length(simul.x))
        uStartDer = zeros(length(simul.y), length(simul.x))


        for i = 1:nIter

            SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

            # L2 norms of the components
            e_0_L2[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uFiltered, u_re, 2)
            e_1_L2[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uDerFiltered, omega * u_im, 2)

            
            grad_e_0_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_e_1_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            grad_c_e_0_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_c_e_1_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            uStart .= simul.uFiltered
            uStartDer .= simul.uDerFiltered 

            println("iter $i of $nIter, $Kx by $Ky elements, exact timestepping")
            
        end



        ########## again with inexact timestepping ##########
        simul.exactStep = false

        uStart = zeros(length(simul.y), length(simul.x))
        uStartDer = zeros(length(simul.y), length(simul.x))

        for i = 1:nIter

            SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

            # L2 norms of the components
            e_0_L2_inexact[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uFiltered, u_re, 2)
            e_1_L2_inexact[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uDerFiltered, omega * u_im, 2)

            grad_e_0_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_e_1_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            grad_c_e_0_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_c_e_1_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            uStart .= simul.uFiltered
            uStartDer .= simul.uDerFiltered 

            println("iter $i of $nIter, $Kx by $Ky elements, inexact timestepping")
            
        end

    end    



    matwrite("WedgeConvHistory.mat", Dict(
    "e_0_L2" => e_0_L2,
    "e_1_L2" => e_1_L2,
    "grad_e_0_L2" => grad_e_0_L2,
    "grad_e_1_L2" => grad_e_1_L2,
    "grad_c_e_0_L2" => grad_c_e_0_L2,
    "grad_c_e_1_L2" => grad_c_e_1_L2,
    "e_0_L2_inexact" => e_0_L2_inexact,
    "e_1_L2_inexact" => e_1_L2_inexact,
    "grad_e_0_L2_inexact" => grad_e_0_L2_inexact,
    "grad_e_1_L2_inexact" => grad_e_1_L2_inexact,
    "grad_c_e_0_L2_inexact" => grad_c_e_0_L2_inexact,
    "grad_c_e_1_L2_inexact" => grad_c_e_1_L2_inexact,
    "u_re_coarse" => u_re_coarse,
    "u_im_coarse" => u_im_coarse,
    "x_coarse" => xVals_coarse,
    "y_coarse" => yVals_coarse,
    "f_coarse" => fVals_coarse,
    "u_re_fine" => u_re_fine,
    "u_im_fine" => u_im_fine,
    "x_fine" => xVals_fine,
    "y_fine" => yVals_fine,
    "f_fine" => fVals_fine,
    "u_re_H1_coarse" => u_re_H1_coarse,
    "u_re_H1_fine" => u_re_H1_fine,
    "u_im_H1_coarse" => u_im_H1_coarse,
    "u_im_H1_fine" => u_im_H1_fine,
    "u_seminorm_coarse" => u_seminorm_coarse,
    "u_seminorm_fine" => u_seminorm_fine,
    "c_data_fine" => c_data_fine
    ))



end


function convHistoryWedgeNewNorm()

    ###### setup ######      
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0






    alphas = [0.0; 0.0; 0.0; 0.0]
    
    omega = 10.0

    nIter = 2000

    N = 6   # degree of interpolation polynomials

    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.

    KxBase = 6
    KyBase = 10
    
    # the data we want: second column has the data from the finer spatial discretisation
    e_0_L2 = zeros(nIter, 2)
    e_1_L2 = zeros(nIter, 2)
    e_0_l2 = zeros(nIter, 2)
    e_1_l2 = zeros(nIter, 2)

    grad_e_0_L2 = zeros(nIter, 2)
    grad_e_1_L2 = zeros(nIter, 2)
    

    grad_c_e_0_L2 = zeros(nIter, 2)
    grad_c_e_1_L2 = zeros(nIter, 2)

    c_e_0_l2 = zeros(nIter, 2)
    c_e_1_l2 = zeros(nIter, 2)

    # the same thing but with inexact timestepping
    e_0_L2_inexact = zeros(nIter, 2)
    e_1_L2_inexact = zeros(nIter, 2)
    e_0_l2_inexact = zeros(nIter, 2)
    e_1_l2_inexact = zeros(nIter, 2)
    
    grad_e_0_L2_inexact = zeros(nIter, 2)
    grad_e_1_L2_inexact = zeros(nIter, 2)
    
    grad_c_e_0_L2_inexact = zeros(nIter, 2)
    grad_c_e_1_L2_inexact = zeros(nIter, 2)
    
    c_e_0_l2_inexact = zeros(nIter, 2)
    c_e_1_l2_inexact = zeros(nIter, 2)



    number_of_refinements = 2

    first_nx = (N-1)*KxBase + 1
    first_ny = (N-1)*KyBase + 1 

    u_re_coarse = zeros(first_nx, first_ny)
    u_im_coarse = zeros(first_nx, first_ny)
    xVals_coarse = zeros(first_nx, 1)
    yVals_coarse = zeros(first_ny, 1)
    fVals_coarse = zeros(first_nx, first_ny)

    second_nx = (N-1)*2*KxBase + 1
    second_ny = (N-1)*2*KyBase + 1

    u_re_fine = zeros(second_nx, second_ny)
    u_im_fine = zeros(second_nx, second_ny)
    xVals_fine = zeros(second_nx, 1)
    yVals_fine = zeros(second_ny, 1)
    fVals_fine = zeros(second_nx, second_ny)
    c_data_fine = zeros(second_nx, second_ny)
    
    u_re_H1_coarse = 0.0
    u_re_H1_fine = 0.0
    u_im_H1_coarse = 0.0
    u_im_H1_fine = 0.0
    u_seminorm_coarse = 0.0
    u_seminorm_fine = 0.0

    L_coarse = zeros(first_nx, first_ny)
    L_fine = zeros(second_nx, second_ny)
    
    M_coarse = zeros(first_nx, first_ny)
    M_fine = zeros(second_nx, second_ny)

    H_coarse = zeros(first_nx, first_ny)
    H_fine = zeros(second_nx, second_ny)

    F_coarse = zeros(first_nx, first_ny)
    F_fine = zeros(second_nx, second_ny)
    
    final_iterate_re = zeros(second_nx, second_ny)
    final_iterate_im = zeros(second_nx, second_ny)

    for j = 1:number_of_refinements

        # choose number of elements
        Kx = Integer(KxBase*j)
        Ky = Integer(KyBase*j)

        # create simul object and forcing
        simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

        nx = length(simul.x)
        ny = length(simul.y)
        DoFs = Int(nx*ny)

        
        #wd = 0.2
        #fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'
        #fVals = 10*exp.(-((simul.y .- 1)/wd).^2) * exp.(-((simul.x .- 0.5)/wd).^2)'
        
        #fVals = 10*pi^2 * cos.(3*pi*simul.y) * cos.(pi*simul.x)'

        fVals = 100 * exp.(-5*((simul.y.-123)./1000).^2) * exp.(-5*((simul.x .- 245.0)./600).^2)'

        g = zeros(length(simul.y), length(simul.x))       


        # compute discretisation of the Laplacian
        L, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, 0.0, alphas, fVals)
        
        M_inv = spdiagm(reshape(1.0 ./ simul.M, DoFs, 1)[:])
        L = -(M_inv*L)
        
        M = spdiagm(reshape(simul.M, DoFs, 1)[:])


        # reference solution
        H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)
        uSolVec = H\F
        uSol = reshape(uSolVec, length(simul.y), length(simul.x))
        u_re = real(uSol)
        u_im = imag(uSol)


        # store the reference solutions for later
        if (j==1)

            u_re_coarse = u_re
            u_im_coarse = u_im

            xVals_coarse = simul.x
            yVals_coarse = simul.y
            fVals_coarse = fVals

            u_re_H1_coarse = SEM_Wave_2d_Updated.SobolevNorm(simul, u_re, 1.0)
            u_im_H1_coarse = SEM_Wave_2d_Updated.SobolevNorm(simul, u_im, 1.0)
            u_seminorm_coarse = SEM_Wave_2d_Updated.SeminormWH(simul, u_re, u_im)
            L_coarse = L

            H_coarse = H
            F_coarse = F
            M_coarse = M
    
        else
            u_re_fine = u_re
            u_im_fine = u_im

            xVals_fine = simul.x
            yVals_fine = simul.y
            fVals_fine = fVals

            u_re_H1_fine = SEM_Wave_2d_Updated.SobolevNorm(simul, u_re, 1.0)
            u_im_H1_fine = SEM_Wave_2d_Updated.SobolevNorm(simul, u_im, 1.0)
            u_seminorm_fine = SEM_Wave_2d_Updated.SeminormWH(simul, u_re, u_im)

            c_data_fine = simul.c_square
            L_fine = L

            H_fine = H
            F_fine = F
            M_fine = M

        end        

        # Waveholtz-appropriate parameters for the wave solver
        Tend = 2*pi/omega
        delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
        delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
        cMax = sqrt(maximum(simul.c_square))
        nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))

        #nsteps = 2*nsteps

        println("using $nsteps timesteps!")

        uStart = zeros(length(simul.y), length(simul.x))
        uStartDer = zeros(length(simul.y), length(simul.x))


        for i = 1:nIter

            SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

            # l2 norms of the components
            e_0_reshaped = reshape(simul.uFiltered - u_re, DoFs, 1)
            e_1_reshaped = reshape(simul.uDerFiltered - omega*u_im, DoFs, 1)


            e_0_l2[i, j] = sqrt((transpose(e_0_reshaped) * e_0_reshaped)[1, 1]) # fulknep 
            e_1_l2[i, j] = sqrt((transpose(e_1_reshaped) * e_1_reshaped)[1, 1])

            c_e_0_l2[i, j] = sqrt((transpose(e_0_reshaped) * (L * e_0_reshaped))[1, 1]) # fulknep 
            c_e_1_l2[i, j] = sqrt((transpose(e_1_reshaped) * (L * e_1_reshaped))[1, 1])


            # L2-norms of the components
            e_0_L2[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uFiltered, u_re, 2)
            e_1_L2[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uDerFiltered, omega * u_im, 2)

            grad_e_0_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_e_1_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            grad_c_e_0_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_c_e_1_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)


            uStart .= simul.uFiltered
            uStartDer .= simul.uDerFiltered 

            println("iter $i of $nIter, $Kx by $Ky elements, exact timestepping")
            
        end

        final_iterate_re = uStart
        final_iterate_im = uStartDer/omega

        ########## again with inexact timestepping ##########
        simul.exactStep = false

        uStart = zeros(length(simul.y), length(simul.x))
        uStartDer = zeros(length(simul.y), length(simul.x))

        for i = 1:nIter

            SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

            # l2 norms of the components
            e_0_reshaped = reshape(simul.uFiltered - u_re, DoFs, 1)
            e_1_reshaped = reshape(simul.uDerFiltered - omega*u_im, DoFs, 1)

            e_0_l2_inexact[i, j] = sqrt((transpose(e_0_reshaped) * e_0_reshaped)[1, 1])
            e_1_l2_inexact[i, j] = sqrt((transpose(e_0_reshaped) * e_0_reshaped)[1, 1])

            c_e_0_l2_inexact[i, j] = sqrt((transpose(e_0_reshaped) * L * e_0_reshaped)[1, 1])
            c_e_1_l2_inexact[i, j] = sqrt((transpose(e_1_reshaped) * L * e_1_reshaped)[1, 1])

            # L2 norms of the components
            e_0_L2_inexact[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uFiltered, u_re, 2)
            e_1_L2_inexact[i, j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uDerFiltered, omega * u_im, 2)

            grad_e_0_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_e_1_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            grad_c_e_0_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_c_e_1_L2_inexact[i, j] = (SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            uStart .= simul.uFiltered
            uStartDer .= simul.uDerFiltered

            println("iter $i of $nIter, $Kx by $Ky elements, inexact timestepping")
            
        end

    end    



    matwrite("WedgeConvHistoryNewNorm.mat", Dict(
    "e_0_L2" => e_0_L2,
    "e_1_L2" => e_1_L2,
    "e_0_l2" => e_0_l2,
    "e_1_l2" => e_1_l2,
    "grad_e_0_L2" => grad_e_0_L2,
    "grad_e_1_L2" => grad_e_1_L2,
    "grad_c_e_0_L2" => grad_c_e_0_L2,
    "grad_c_e_1_L2" => grad_c_e_1_L2,
    "c_e_0_l2" => c_e_0_l2,
    "c_e_1_l2" => c_e_1_l2,
    "e_0_L2_inexact" => e_0_L2_inexact,
    "e_1_L2_inexact" => e_1_L2_inexact,
    "e_0_l2_inexact" => e_0_l2_inexact,
    "e_1_l2_inexact" => e_1_l2_inexact,
    "grad_e_0_L2_inexact" => grad_e_0_L2_inexact,
    "grad_e_1_L2_inexact" => grad_e_1_L2_inexact,
    "grad_c_e_0_L2_inexact" => grad_c_e_0_L2_inexact,
    "grad_c_e_1_L2_inexact" => grad_c_e_1_L2_inexact,
    "c_e_0_l2_inexact" => c_e_0_l2_inexact,
    "c_e_1_l2_inexact" => c_e_1_l2_inexact,
    "u_re_coarse" => u_re_coarse,
    "u_im_coarse" => u_im_coarse,
    "x_coarse" => xVals_coarse,
    "y_coarse" => yVals_coarse,
    "f_coarse" => fVals_coarse,
    "u_re_fine" => u_re_fine,
    "u_im_fine" => u_im_fine,
    "u_re_whi" => final_iterate_re,
    "u_im_whi" => final_iterate_im,
    "x_fine" => xVals_fine,
    "y_fine" => yVals_fine,
    "f_fine" => fVals_fine,
    "u_re_H1_coarse" => u_re_H1_coarse,
    "u_re_H1_fine" => u_re_H1_fine,
    "u_im_H1_coarse" => u_im_H1_coarse,
    "u_im_H1_fine" => u_im_H1_fine,
    "u_seminorm_coarse" => u_seminorm_coarse,
    "u_seminorm_fine" => u_seminorm_fine,
    "c_data_fine" => c_data_fine,
    "L_coarse" => L_coarse,
    "L_fine" => L_fine,
    "H_coarse" => H_coarse,
    "H_fine" => H_fine,
    "F_coarse" => F_coarse,
    "F_fine" => F_fine,
    "M_coarse" => M_coarse,
    "M_fine" => M_fine
    ))

end

function WedgePlotData()

    ###### setup ######      
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0


    alphas = [0.0; 0.0; 0.0; 0.0]
    
    
    omega = 100.0
    tol = 1e-9

    N = 10   # degree of interpolation polynomials

    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.

    Kx = 18
    Ky = 30


    # create simul object and forcing
    simul_c = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    
    g = zeros(length(simul_c.y), length(simul_c.x))       
    fVals_coarse = 100 * exp.(-5*((simul_c.y.-123)./1000).^2) * exp.(-5*((simul_c.x .- 245.0)./600).^2)'

    println("coarse run")
    u_re_coarse, u_im_coarse, hist = SEM_Wave_2d_Updated.WaveholtzGMRES(simul_c, omega, fVals_coarse, alphas, g, tol)
    u_im_coarse = u_im_coarse/omega

    # store some stuff for later
    c_data_fine = simul_c.c_square    
    xVals_coarse = simul_c.x
    yVals_coarse = simul_c.y

    # create simul object and forcing, finer run
    simul_f = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], 2*Kx, 2*Ky, N, c_square)

    # data for WaveHoltz
    g = zeros(length(simul_f.y), length(simul_f.x))       
    fVals_fine = 100 * exp.(-5*((simul_f.y.-123)./1000).^2) * exp.(-5*((simul_f.x .- 245.0)./600).^2)'

    println("fine run")
    println("fine run not carried out")
    #u_re_fine, u_im_fine, hist = SEM_Wave_2d_Updated.WaveholtzGMRES(simul_f, omega, fVals_fine, alphas, g, tol)
    #u_im_fine = u_im_fine/omega
    u_re_fine = zeros(length(simul_f.y), length(simul_f.x))
    u_im_fine = zeros(length(simul_f.y), length(simul_f.x))

    # store some stuff for later
    xVals_fine = simul_f.x
    yVals_fine = simul_f.y
    c_data_fine = simul_f.c_square    


    matwrite("WedgePlotData.mat", Dict(
    "u_re_coarse" => u_re_coarse,
    "u_im_coarse" => u_im_coarse,
    "x_coarse" => xVals_coarse,
    "y_coarse" => yVals_coarse,
    "f_coarse" => fVals_coarse,
    "c_data_fine" => c_data_fine,
    "u_re_fine" => u_re_fine,
    "u_im_fine" => u_im_fine,
    "x_fine" => xVals_fine,
    "y_fine" => yVals_fine,
    "f_fine" => fVals_fine,
    "c_data_fine" => c_data_fine,
    "omega" => omega
    ))

end