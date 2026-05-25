

using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using Base.Threads

include("SEM_Wave_2d_Updated.jl")
using .SEM_Wave_2d_Updated

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
    Kx = 8 # number of elements in x-direction
    Ky = 8 # number of elements in y-direction

    c_square(x, y) = 1

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    #fVals = zeros(length(simul.y), length(simul.x))
    fVals = 1e3*exp.(-((simul.y .- 0.2*(yr+yl))/0.03).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.03).^2)'

    #alphas = [0; 1/sqrt(2); 0; 1/sqrt(2)] #impedance up and down, Neumann left and right
    alphas = [0; 0.0; 0; 1/sqrt(2)] #impedance down, else Neumann 
    #alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] #impedance all around
    #alphas = [0.0; 0.0; 0.0; 0.0] #Neumann all around
    

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    omega = 14.1

    Tend = 5.0
    nsteps = 1000

    # construct initial data
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    #uStart = exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'
    #uStart = 2 .- (abs.(simul.y) .+ abs.(simul.x'))

    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    #uStartDer = 1e3*exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'

    SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, true, 1, 1.5)


end

function WaveholtzTest()

        
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0


    N = 6   # degree of interpolation polynomials
    Kx = 12 # number of elements in x-direction
    Ky = 12 # number of elements in y-direction

    c_square(x, y) = 1 - 0.999 .* exp.(-((y .- 0.7)/0.1).^2) * exp.(-((x .- 0.3)/0.1).^2)'

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    #fVals = zeros(length(simul.y), length(simul.x))
    fVals = 1e3*exp.(-((simul.y .- 0.2*(yr+yl))/0.03).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.03).^2)'

    #alphas = [0; 1/sqrt(2); 0; 1/sqrt(2)] #impedance up and down, Neumann left and right
    #alphas = [0; 0; 0; 1/sqrt(2)] #impedance down, else Neumann 
    #alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] #impedance all around
    alphas = [0.0; 0.0; 0.0; 0.0] #Neumann all around
    #alphas = [1/sqrt(2); 1/sqrt(2); 0; 1/sqrt(2)] #Neumann left, else impedance
    #alphas = [1/sqrt(2); 0; 0; 1/sqrt(2)] #Neumann left and up, else impedance
    

    g = zeros(length(simul.y), length(simul.x))

    simul.g = g
    omega = 15.123


    # construct initial data
    uStart = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    #uStart = exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'
    #uStart = 2 .- (abs.(simul.y) .+ abs.(simul.x'))

    uStartDer = zeros(length(simul.y), 1) * zeros(1, length(simul.x))
    #uStartDer = 1e3*exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.1*(xr+xl))/0.1).^2)'

    #uSol, nIter = SEM_Wave_2d_Updated.Waveholtz(simul, omega, fVals, alphas, g, 1e-2)


    #=
    uSol, u_1, historyExact = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, 1e-3)

    plt = heatmap(simul.x, simul.y, uSol)
    savefig(plt, "WaveholtzTestExact")

    gmresIterExact = length(historyExact.data[:resnorm])
    

    simul.exactStep = false

    uSol, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, 1e-3)

    plt = heatmap(simul.x, simul.y, uSol)
    savefig(plt, "WaveholtzTest")


    plt = plot(LinRange(1, gmresIterExact, gmresIterExact), historyExact.data[:resnorm], yscale=:log10, label = "WaveHoltz with Exact Scheme")

    gmresIter = length(history.data[:resnorm])
    plot!(LinRange(1, gmresIter, gmresIter), history.data[:resnorm], yscale=:log10, label = "Normal WaveHoltz")
    savefig(plt, "WaveholtzConvComp")

    =#

    maxIter = 100
    u_0, dataExact = SEM_Wave_2d_Updated.WaveholtzConvHistory(simul, omega, fVals, alphas, g, maxIter)

    simul.exactStep = false
    u_0, data = SEM_Wave_2d_Updated.WaveholtzConvHistory(simul, omega, fVals, alphas, g, maxIter)

    plt = plot(LinRange(1, maxIter, maxIter), dataExact, xscale=:log10, yscale=:log10, label = "WaveHoltz with Exact Scheme")
    plot!(LinRange(1, maxIter, maxIter), data, xscale=:log10, yscale=:log10, label = "Normal WaveHoltz iteration")
    savefig(plt, "WHIConvComp")

end



function frequencySweepBox()
    
    ###### setup ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    #alphas = [0.0; 0.0; 0.0; 1/sqrt(2)] #impedance down, else Neumann 
    alphas = [0; 1/sqrt(2); 0; 1/sqrt(2)] #impedance up and down, Neumann left and right
    #alphas = [1/sqrt(2); 1/sqrt(2); 0; 1/sqrt(2)] #Neumann left, else impedance
    #alphas = [1/sqrt(2); 0; 0; 1/sqrt(2)] #Neumann left and up, else impedance

    spaceTol = 1e-3
    systemTol = 1e-6

    
    h = 0.05
    omegaSmall = 2.0
    omegaBig = 15.0
    
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
                    
                    #println("thread $tid considering value $omega")

                    current_omega[tid] = omega

                    N = 6   # degree of interpolation polynomials

                    c_square(x, y) = 1

                    L = maximum([xr-xl, yr-yl]) # length of the domain

                    lambda = 1/(omega/(2*pi))
                    N_lambda = L / lambda
                    
                    # rule of thumb 
                    PPW = pi*(N_lambda/spaceTol)^(1/N)

                    dofs_x = ceil(PPW * N_lambda)       # how many points do we want per dimension? points = point/wavelength * wavelength

                    Kx = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements
                    Ky = Kx # the same number in both directions!

                     
                    
                    # if one wants a floor on the number of elements
                    
                    if Kx < 3
                        Kx = 3
                        Ky = Kx
                    end
                    

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
                    fVals = 10*exp.(-((simul.y .- 0.223)/0.01).^2) * exp.(-((simul.x .- 0.652)/0.01).^2)'
                    g = zeros(length(simul.y), length(simul.x))       

                    
                    

                    println("$omega: $Kx by $Ky")
                    
                    sol1, nIter = SEM_Wave_2d_Updated.Waveholtz(simul, omega, fVals, alphas, g, systemTol)

                    
                    

                    #println("gmres || || omega: " * string(omega) * " || " * "number of elements: " * string(Ky) * " by " * string(Kx))
                    #println("thread $tid started gmres solve for omega = $omega")

                    u_0, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas, g, systemTol)

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

    savefig(plt, "BoxIterPlot")

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




function mmsTestQuick()
    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    #=
    xl = -2.123/5
    xr = 4.0234/5
    yl = -3.234567/5
    yr = 1.5678/5
    =#

    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0




    N = 8   # degree of interpolation polynomials
    Kx = 16 # number of elements in x-direction
    Ky = 16 # number of elements in y-direction

    

    c_square(x, y) = 1
    omega = 7.1
    Tend = 1.0
    nsteps = 500

 
    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    simul.useMMS = true # we mean to perform an mms test


    #=
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
    =#


    # trigonometric test
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
    



    ###### set up the problem ######
    fVals = 1e4*exp.(-((simul.y .- 0.2*(yr+yl))/0.03).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.03).^2)' # should not matter what we put here, as the MMS functionality should put the proper values for f regardless



    alphas = [1/sqrt(2); 1/sqrt(2); 0; 1/sqrt(2)] #Neumann left, else impedance

    g = zeros(length(simul.y), length(simul.x))


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

    SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, true, 10, 5)

    plt = surface(simul.x, simul.y[end:-1:1], log10.(abs.(simul.uNow-reference)))
    savefig(plt, "log10diff")

end


function mmsConvPlot()
    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######
    
    
    xl = -2.123/5
    xr = 4.0234/5
    yl = -3.234567/5
    yr = 1.5678/5
    

    
    #xl = -1.0
    #xr = 1.0
    #yl = -1.0
    #yr = 1.0
    

    N = 3 # degree of interpolation polynomials
    numberOfRuns = 5
    errorData = zeros(numberOfRuns, 1)
    KxVec = zeros(numberOfRuns, 1)
    
    alphas = [0.0; 0.0; 0.0; 0.0] #Neumann everywhere
    alphas = [1/sqrt(2); 1/sqrt(2); 0; 1/sqrt(2)] #Neumann left, else impedance
    #alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] #impedance everywhere
    alphas = rand(4) #random alphas

    for j = numberOfRuns:-1:1

        c_square(x, y) = 1
        omega = 7.1
        Tend = 0.0234
        nsteps = Int(5*2^j)

        Kx = Int(2^(j)) # number of elements in x-direction
        Ky = Kx # number of elements in y-direction

        KxVec[j] = Kx

        println("$Kx by $Ky elements running")

    
        simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

        simul.useMMS = true # we mean to perform an mms test


        #=
        # polynomial test
        simul.MMS_j.type = [3 3 3]

        simul.MMS_j.coeff[1, 1] = 1
        #simul.MMS_j.coeff[1, 2] = 0.567
        simul.MMS_j.coeff[1, 2] = 0.0
        simul.MMS_j.coeff[1, 3] = 1.0

        simul.MMS_j.coeff[2, 1] = 1
        simul.MMS_j.coeff[2, 2] = 0.0
        simul.MMS_j.coeff[2, 3] = 0.0

        simul.MMS_j.coeff[3, 1] = 1
        simul.MMS_j.coeff[3, 2] = 0
        simul.MMS_j.coeff[3, 3] = 4
        =#


        
        # trigonometric test
        simul.MMS_j.type = [2 2 2]
        simul.MMS_j.coeff[1, 1] = 1.0
        simul.MMS_j.coeff[1, 2] = 2*pi
        simul.MMS_j.coeff[1, 3] = 0.0
        
        simul.MMS_j.coeff[2, 1] = 1
        simul.MMS_j.coeff[2, 2] = 2*pi
        simul.MMS_j.coeff[2, 3] = 0.0

        simul.MMS_j.coeff[3, 1] = 1.0
        simul.MMS_j.coeff[3, 2] = 2.0
        simul.MMS_j.coeff[3, 3] = 0.0



        # trigonometric test
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
        
        


        ###### set up the problem ######
        fVals = 1e4*exp.(-((simul.y .- 0.2*(yr+yl))/0.03).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.03).^2)' # should not matter what we put here, as the MMS functionality should put the proper values for f regardless


        g = zeros(length(simul.y), length(simul.x))


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

        SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)

        errorData[j] = maximum(abs.(simul.uNow - reference))

        #plt = surface(simul.x, simul.y[end:-1:1], abs.(simul.uNow - reference))    
        plt = surface(simul.x, simul.y[end:-1:1], simul.uNow - reference)    
        #plt = surface(simul.x, simul.y[end:-1:1], reference)
        savefig(plt, "testsurface" * string(j))

        println(maximum(reference))
        println(maximum(simul.uNow))


    end

    println(errorData)
    plt = plot(KxVec, errorData, xscale=:log10, yscale=:log10)
    savefig(plt, "mmsConvPlot")

end


function biggerMMSTest()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######

    
    xl = -2.123/5
    xr = 4.0234/5
    yl = -3.234567/5
    yr = 1.5678/5
    

    numberOfRuns = 5
    numbersOfElements = zeros(numberOfRuns)

    for i = 1:numberOfRuns

        #numbersOfElements[i] = Integer(2^(i+2))
        numbersOfElements[i] = Integer(round(2^((i+1)/3+3)))
        
    end


    Ns = [2, 4, 6]

    plt = plot()


    alphas = rand(4) #random alphas

    for i = length(Ns):-1:1

        N = Ns[i]  # degree of interpolation polynomials
        
        c_square(x, y) = 1

        diffs = zeros(numberOfRuns)

        for j = numberOfRuns:-1:1

            println(numbersOfElements[j])
            println(N)

            c_square(x, y) = 1
            omega = 7.1
            Tend = 0.00234
            nsteps = Integer(10*numbersOfElements[j])
            #nsteps = 10
            println("steps: $nsteps")

            Kx = Integer(numbersOfElements[j]) # number of elements in x-direction
            Ky = Kx # number of elements in y-direction

            simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

            simul.useMMS = true # we mean to perform an mms test

            fVals = zeros(length(simul.y), length(simul.x))
            g = zeros(length(simul.y), length(simul.x))


            uStart = exp.(-((simul.y .- 0.5*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.5*(xr+xl))/0.1).^2)'
            uStartDer = zeros(length(simul.y), length(simul.x))

            reference = zeros(length(simul.y), length(simul.x))




            # trigonometric test
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


            #=
            # polynomial test
            simul.MMS_j.type = [3 3 3]

            simul.MMS_j.coeff[1, 1] = 1
            simul.MMS_j.coeff[1, 2] = 0.567
            #simul.MMS_j.coeff[1, 2] = 0.0
            simul.MMS_j.coeff[1, 3] = 2.0

            simul.MMS_j.coeff[2, 1] = sqrt(11)
            simul.MMS_j.coeff[2, 2] = (sqrt(5) + 1)/2
            simul.MMS_j.coeff[2, 3] = 5.0

            simul.MMS_j.coeff[3, 1] = 1
            simul.MMS_j.coeff[3, 2] = 0
            simul.MMS_j.coeff[3, 3] = 1
            =#


            if simul.useMMS
                for s = 1:length(simul.y)
                    for t = 1:length(simul.x)
                        uStart[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 0, simul.MMS_j)
                        uStartDer[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], 0.0, 0, 0, 1, simul.MMS_j)
                        reference[s, t] = MMS.MMSfun(simul.x[t], simul.y[end+1-s], Tend, 0, 0, 0, simul.MMS_j)
                    end
                end
            end


            # decide whether you want the exact timestepping for WaveHoltz
            #simul.exactStep = false

            SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)
            
            diffs[j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uNow,reference, 2) / SEM_Wave_2d_Updated.LpNorm(simul, simul.uNow, zeros(length(simul.y), length(simul.x)), 2)
            #diffs[j] = maximum(abs.(simul.uNow - reference))

        end

        println(diffs)

        logLine = zeros(numberOfRuns, 1)
        
        #=
        logLine[end] = diffs[end]

        for j = numberOfRuns-2:-1:1
            logLine[j] = logLine[j+1] * (numbersOfElements[j+1] / numbersOfElements[j])^(N+1)
        end
        =#


        logLine[1] = diffs[1]

        logLineOrder = N
        #logLineOrder = 1
        for j = 2:numberOfRuns
            logLine[j] = logLine[j-1] * (numbersOfElements[j-1] / numbersOfElements[j])^logLineOrder
        end


        scatter!(numbersOfElements, diffs, xscale=:log10, yscale=:log10, label = "N = " * string(N), ms = 3)
        plot!(numbersOfElements, logLine, xscale=:log10, yscale=:log10, label = "reference, order = $logLineOrder", linestyle=:dash, legend=:bottomleft)
    

    end
    
    savefig(plt, "mmsConvPlotMixed")

end



function selfConvergenceTest()

    
    ###### set up the time interval, the number of elements, the degree of the interpolation polynomials ######

    
    xl = -0.01
    xr = 0.01
    yl = -0.01
    yr = 0.01
    

    #=
    xl = -1.0
    xr = 1.0
    yl = -1.0
    yr = 1.0
    =#

    numberOfRuns = 4
    numbersOfElements = zeros(numberOfRuns)

    for i = 1:numberOfRuns

        numbersOfElements[i] = Integer(2^(i+2))
        #numbersOfElements[i] = Integer(round(2^((i+1)/2+1)))
        
    end


    Ns = [1, 3, 5]

    plt = plot()

    for i = length(Ns):-1:1

        N = Ns[i]  # degree of interpolation polynomials
        alphas = [0.0; 0.0; 0.0; 0.0] #Neumann everywhere
        
        #c_square(x, y) = 1

        
        c_square(x, y) = 1 + (x-0.2).^2


        reference = zeros(Integer(numbersOfElements[end] * N + 1), (Integer(numbersOfElements[end] * N + 1)))

        diffs = zeros(numberOfRuns-1)

        for j = numberOfRuns:-1:1

            println(numbersOfElements[j])


            Kx = Integer(numbersOfElements[j]) # number of elements in x-direction
            Ky = Integer(numbersOfElements[j]) # number of elements in y-direction

            simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)


            #fVals = exp.(-((simul.y .- 0.2*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.7*(xr+xl))/0.1).^2)'
            fVals = zeros(length(simul.y), length(simul.x))
            g = zeros(length(simul.y), length(simul.x))

            omega = 7.1

            Tend = 0.00052345
            nsteps = 50

            uStart = exp.(-((simul.y .- 0.5*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.5*(xr+xl))/0.1).^2)'

            uStartDer = zeros(length(simul.y), length(simul.x))

            SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)


            ### Compare to finest run ###

            if (j == length(numbersOfElements))
                reference = simul.uNow
            else    

                coarserReference = reference[1:N*2^(length(numbersOfElements) - j):end, 1:N*2^(length(numbersOfElements) - j):end]

                diffs[j] = maximum(abs.(simul.uNow[1:N:end, 1:N:end] - coarserReference)./abs.(coarserReference))
                #diffs[j] = maximum(abs.(reference[1:N*2^(length(numbersOfElements) - j):end, 1:N*2^(length(numbersOfElements) - j):end] - simul.uNow[1:N:end, 1:N:end]))
                
            end

        end

        println(diffs)

        logLine = zeros(numberOfRuns-1, 1)
        
        #=
        logLine[end] = diffs[end]

        for j = numberOfRuns-2:-1:1
            logLine[j] = logLine[j+1] * (numbersOfElements[j+1] / numbersOfElements[j])^(N+1)
        end
        =#


        logLine[1] = diffs[1]

        logLineOrder = N+1
        #logLineOrder = 1
        for j = 2:numberOfRuns-1
            logLine[j] = logLine[j-1] * (numbersOfElements[j-1] / numbersOfElements[j])^logLineOrder
        end


        println(numbersOfElements[1:end-1])
        scatter!(numbersOfElements[1:end-1], diffs, xscale=:log10, yscale=:log10, label = "N = " * string(N), ms = 3)
        plot!(numbersOfElements[1:end-1], logLine, xscale=:log10, yscale=:log10, label = "reference, order = $logLineOrder", linestyle=:dash, legend=:bottomleft)
    

    end
    
    savefig(plt, "selfConvPlot")

end



function neumannComparison()
    # checks whether the new and old SEM_Wave-modules agree for nonmixed BCs 

    xl = -2.123/5
    xr = 4.0234/5
    yl = -3.234567/5
    yr = 1.5678/5


    #a = 1/sqrt(2)
    a = 0
    b = sqrt(1-a^2)
    bc = [a, b]

    omega = 7.1
    Tend = 0.0234
    nsteps = 25

    Kx = 10
    Ky = Kx
    c_square(x, y) = 1

    N = 3
    
    simul = SEM_Wave_2d.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
    simulNew = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
    
    
    fVals = 10*exp.(-((simul.y .- 0.223)/0.01).^2) * exp.(-((simul.x .- 0.652)/0.01).^2)'
    g = zeros(length(simul.y), length(simul.x))       




    uStart = exp.(-((simul.y .- 0.5*(yr+yl))/0.1).^2) * exp.(-((simul.x .- 0.5*(xr+xl))/0.1).^2)'
    uStartDer = exp.(-((simul.y .- 0.5*(yr+yl))/0.2).^2) * exp.(-((simul.x .- 0.2*(xr+xl))/0.1).^2)'


    alpha = 0.1
    bc = [alpha, sqrt(1-alpha^2)]
    alphas = [alpha; alpha; alpha; alpha]

    
    #=
    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, false)
    SEM_Wave_2d_Updated.Simulate(simulNew, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)

    
    println(maximum(abs.(simul.uNow - simulNew.uNow)))
    plt = surface(simul.x, simul.y[end:-1:1], simul.uNow - simulNew.uNow)    
    println((simul.uNow - simulNew.uNow)[end-5:end, end-5:end])
    =#

    #=
    # initialisation test
    SEM_Wave_2d_Updated.Initialise!(simulNew, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g)
    SEM_Wave_2d.Initialise!(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g)
    testing = simul.uPrev
    testingNew = simulNew.uPrev
    =#

    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, false)
    SEM_Wave_2d_Updated.Simulate(simulNew, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)
    testing = simul.uNow
    testingNew = simulNew.uNow

    #simulNew.uNow = simul.uNow
    #simulNew.uPrev = simul.uPrev


    
    #testing = SEM_Wave_2d.TimeSteppingTerm(simul)
    #testingNew = SEM_Wave_2d_Updated.TimeSteppingTerm(simulNew)

    #testing = SEM_Wave_2d.LaplaceTerm(simul, simul.uNow)
    #testingNew = SEM_Wave_2d_Updated.LaplaceTerm(simulNew, simul.uNow)

    #testing = SEM_Wave_2d.ForcingTerm(simul, nsteps + 1)
    #testingNew = SEM_Wave_2d_Updated.ForcingTerm(simulNew, nsteps + 1)
    

    println(maximum(abs.(testing - testingNew)))

    plt = surface(simul.x, simul.y[end:-1:1], testing - testingNew)
    savefig(plt, "neumannComparison")




end
