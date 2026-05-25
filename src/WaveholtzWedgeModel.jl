


using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using Base.Threads
using Dates

include("SEM_Wave_2d.jl")
using .SEM_Wave_2d


include("MMS.jl")
using .MMS




function main()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    

    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0

    
    ###### set up the problem ######
    piMultiples = [10.0, 20.0, 40.0]
    a = 1/sqrt(2)
    b = sqrt(1-a^2)
    bc = [a, b]
    tol = 1e-7

    iterCount = 20 # number of iterations we perform for our convergence history

    # aja baja!
    #=
    delta(x) = Float64(x==0.0) # the best implementation known to man
    f(x, y) = omega^2 * delta(abs(x - 300.0))' * delta(abs(y))
    fVals = f.(simul.x', simul.y)
    =#

    # step 1: find values for (r*, s*) corresponding to the location of the source (x*, y*) 
    # step 2: evaluate the basis functions at (r*, s*)
    # step 3: integral over r, (-1,1) times integral over s (-1, 1), integrand f(x(r, s), y(r, s))*\psi(r, s) drds
    #  => omega^2 [Jacobian] \psi(r*, s*)

    finestXvals = NaN
    finestYvals = NaN
    uReal = NaN

    figure = plot(legend=:bottomleft)

    for j = 1:length(piMultiples)

        omega = piMultiples[j] * pi

        N = 4   # degree of interpolation polynomials
        c_min = 1000 # has to be put together
        L = maximum([xr-xl, yr-yl]) # max side length of the domain

        # rule of thumb for 4th order, rounded up a bit
        dofs_x = Int(ceil(1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
        Kx = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!

        Ky = Integer(ceil(Kx * (5/3))) # number of elements in y-direction

        println("frequency: " * string(piMultiples[j]) * "π" * " || elements: " * string(Ky) * " by " * string(Kx))
        

        c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.

        simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

        fVals = omega^2 * exp.(-0.1*(simul.y).^2) * exp.(-0.1*(simul.x .- 300.0).^2)'

        g = zeros(length(simul.y), length(simul.x))

        u_WHI, data = SEM_Wave_2d.WaveholtzConvHistory(simul, omega, fVals, bc, g, iterCount, false)
        
        plot!(data, yscale=:log10, label = "errors at frequency " * string(piMultiples[j]) * "π")

        if j == length(piMultiples)
            finestXvals = simul.x
            finestYvals = simul.y
            uReal = u_WHI[1:length(simul.y), 1:end]
        end

    end

    savefig(figure, "WedgeModelConvTest")
    
    plt = heatmap(finestXvals, finestYvals, log10.(abs.(uReal)), aspect_ratio =:1.00)
    savefig(plt, "WaveholtzWedge")

    #u_0, u_1, history = SEM_Wave_2d.WaveholtzGMRESnew(simul, omega, fVals, bc, g, tol)    
    #println(history.data[:resnorm])

    #plt = plot(history.data[:resnorm], yscale=:log10)
    #savefig(plt, "gmresConvHistoryWedge")
    
    #plt = heatmap(simul.x, simul.y, log10.(abs.(u_WHI)), aspect_ratio =:1.00)
    
    
end


function frequencySweep()

    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0


    #a = 1/sqrt(2)
    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    
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
    omegaSmall = 3.9 + h
    omegaBig = 4.95
    
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


                    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
                    
                    fVals = 100 * exp.(-5*((simul.y.-123)./1000).^2) * exp.(-5*((simul.x .- 245.0)./600).^2)'

                    
                    g = zeros(length(simul.y), length(simul.x))        

                    
                    

                    println("$omega: $Kx by $Ky")
                    
                    sol1, nIter = SEM_Wave_2d.Waveholtz(simul, omega, fVals, bc, g, systemTol)

                    
                    

                    #println("gmres || || omega: " * string(omega) * " || " * "number of elements: " * string(Ky) * " by " * string(Kx))
                    #println("thread $tid started gmres solve for omega = $omega")

                    u_0, u_1, history = SEM_Wave_2d.WaveholtzGMRESnew(simul, omega, fVals, bc, g, systemTol)

                    gmresIter = length(history.data[:resnorm])
                    #println("thread $tid finished gmres solve for omega = $omega in $gmresIter iterations")

                    
                    gmresIters[j] = gmresIter
                    whIters[j] = nIter


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



function makeAPlot()

    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    xl = 0.0
    xr = 600.0
    yl = 0.0
    yr = 1000.0


    #a = 1/sqrt(2)
    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    spaceTol = 1e-7
    systemTol = 1e-7

    #omega = 15*pi
    #omega = 66.0 # about two wavelengths in the first part
    #omega = 5.0

    omega = 4.3 # close to the first nonzero eigenmode
    omega = 9.275 # close to the second nonzero eigenmode
    omega = 11.35 # close to the third nonzero eigenmode
    omega = 12.6 # close to the fourth nonzero eigenmode

    N = 6   # degree of interpolation polynomials
    
    # variable wave speed for the problem
    c_square(x, y) = (2100 - 1100*Float64(y > (x/6 + 400)) + 1900*Float64(y > (800 - x/3))).^2 # wedge model.
    c_min = 1000 # has to be put correctly

    ##################################### Copied from frequencySweep() #####################################

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


    println("omega: $omega with $Kx by $Ky elements")
    

    #=          older choice of elements
    L = maximum([xr-xl, yr-yl]) # max side length of the domain

    # rule of thumb for 4th order, rounded up a bit
    dofs_x = Int(ceil(2.5*1.8 * (L * omega / (c_min * tol))^(1/N))) # how many points do we want per dimension?
    Kx = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements!

    Ky = Integer(ceil(Kx * (5/3))) # number of elements in y-direction
    
    =#

    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
    fVals = omega^2 * exp.(-0.1*(simul.y).^2) * exp.(-0.1*(simul.x .- 300.0).^2)'

    
    g = zeros(length(simul.y), length(simul.x))        


    #println("WHI || || omega: " * string(omega) * " || " * "number of elements: " * string(Ky) * " by " * string(Kx))
    #sol, nIter = SEM_Wave_2d.Waveholtz(simul, omega, fVals, bc, g, tol)
    

    #sol, data = SEM_Wave_2d.WaveholtzConvHistory(simul, omega, fVals, bc, g, 100, false)
    #plt = plot(data, yscale=:log10, label = "errors at frequency " * string(omega)) 
    #savefig(plt, "wedgeStrangeConvHist")


    u_0, u_1, history = SEM_Wave_2d.WaveholtzGMRESnew(simul, omega, fVals, bc, g, systemTol)
    gmresIter = length(history.data[:resnorm])

    println("waveholtz with gmres finished in $gmresIter iterations")


    # plot the result
    plt = heatmap(simul.x, simul.y, u_0)
    savefig(plt, "WedgePlot")

end

