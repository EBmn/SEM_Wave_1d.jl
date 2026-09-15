


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


function frequencySweepBox()
    
    ###### setup ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    
    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)]
    #alphas = [0.0; 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)]
    #alphas = [0.0; 1/sqrt(2); 1/sqrt(2); 0.0]
    #alphas = [1/sqrt(2); 0.0; 1/sqrt(2); 0.0]
    #alphas = [1/sqrt(2); 0.0; 0.0; 0.0]
    #alphas = [0.0; 0.0; 0.0; 0.0]


    spaceTol = 1e-3
    systemTol = 1e-6

    
    #h = 0.05
    #omegaSmall = 15.0 + h
    #omegaBig = 20.0
    h = 0.05
    omegaSmall = 2.0
    omegaBig = 15.0

    numberOfOmegas = Int(round((omegaBig-omegaSmall) / h + 1))
    omegas = collect(LinRange(omegaSmall, omegaBig, numberOfOmegas))


    #omegas = [3.0, 3.012872059417254, 3.025744118834508, 3.0386161782517616, 3.0514882376690156, 3.0643602970862696, 3.0772323565035236, 3.0901044159207776, 3.1029764753380316, 3.115848534755285, 3.128720594172539, 3.259891770368754, 3.3781908871477153, 3.4964900039266764, 3.614789120705638, 3.733088237484599, 3.85138735426356, 3.9696864710425213, 4.087985587821483, 4.206284704600444, 4.324583821379405, 4.6101831535239315, 4.777483368889497, 4.944783584255062, 5.112083799620628, 5.279384014986194, 5.446684230351758, 5.613984445717325, 5.78128466108289, 5.9485848764484555, 6.115885091814021, 6.350606163894235, 6.418027020608885, 6.485447877323534, 6.552868734038182, 6.620289590752832, 6.687710447467481, 6.755131304182131, 6.822552160896779, 6.889973017611428, 6.957393874326078, 7.193992107884, 7.3631694847272735, 7.532346861570547, 7.701524238413819, 7.870701615257093, 8.039878992100366, 8.209056368943639, 8.378233745786913, 8.547411122630185, 8.71658849947346, 8.934766974903336, 8.98376807348994, 9.032769172076545, 9.08177027066315, 9.130771369249754, 9.179772467836358, 9.228773566422962, 9.277774665009566, 9.32677576359617, 9.375776862182775, 9.471124352135446, 9.51747074350151, 9.563817134867577, 9.610163526233642, 9.656509917599708, 9.702856308965773, 9.749202700331839, 9.795549091697904, 9.84189548306397, 9.888241874430035, 10.061186914281818] #, 10.187785562767534, 10.314384211253248, 10.440982859738964, 10.56758150822468, 10.694180156710397, 10.820778805196113, 10.947377453681828, 11.073976102167544, 11.20057475065326, 11.439827691431722, 11.552481983724466, 11.665136276017211, 11.777790568309957, 11.890444860602702, 12.003099152895446, 12.115753445188192, 12.228407737480937, 12.341062029773681, 12.453716322066427, 12.60152949881881, 12.636688383278448, 12.671847267738086, 12.707006152197724, 12.742165036657362, 12.777323921117, 12.812482805576638, 12.847641690036276, 12.882800574495914, 12.917959458955552, 12.987257477147908, 13.05631333650028, 13.12536919585265, 13.19442505520502, 13.26348091455739, 13.41129409130978, 13.48034995066215, 13.54940581001452, 13.61846166936689, 13.68751752871926, 13.82411685087773, 13.9607161762443, 14.09731550161087, 14.23391482697744, 14.37051415234401, 14.50711347771058, 14.64371280307715, 14.78031212844372, 14.91691145381029, 15.05351077917686, 15.19011010454343, 15.2591659639058, 15.32822182326817, 15.39727768263054, 15.46633354199291, 15.53538940135528, 15.60444526071765, 15.77361258418453, 15.8426684435469, 15.91172430290927, 15.98078016227164, 16.083433554221595, 16.18608694617155, 16.288740338121505, 16.39139373007146, 16.494047122021415, 16.59670051397137, 16.69935390592133, 16.802007297871285, 16.90466068982124, 17.07382801106811, 17.17648140301807, 17.27913479496802, 17.38178818691798, 17.48444157886794, 17.5870949708179, 17.68974836276785, 17.79240175471781, 17.89505514666777, 17.99770853861773, 18.08004530056578, 18.162382062513835, 18.24471882446189, 18.32705558640995, 18.47486490856852, 18.557191670516575, 18.63951843246463, 18.721845194412685, 18.80417195636074, 18.8864987183088, 18.94881958056486, 18.97015881081482, 19.004838595331405, 19.039518379790987, 19.07419816425057, 19.10887794871015, 19.247673797420887, 19.316729656773255, 19.385785516125623, 19.45484137547799, 19.52389723483036, 19.59295309418273, 19.662008953535096, 19.731064812887464, 19.800120672239835, 19.891615761842164, 19.914054992092126, 19.936494222342088, 19.95893345259205, 19.98137268284201, 20.003811913091976, 20.026251143341938, 20.0486903735919, 20.07112960384186, 20.093568834091823, 20.203138621503637, 20.290269178665493, 20.377399735827346, 20.4645302929892, 20.551660850151055, 20.638791407312908, 20.725921964474765, 20.813052521636617, 20.90018307879847, 20.987313635960327, 21.15778095512276, 21.241117717123338, 21.32445447912392, 21.407791241124496, 21.491128003125077, 21.574464765125654, 21.657801527126235, 21.741138289126813, 21.824475051127393, 21.90781181312797, 22.011445494734303, 22.03174241434006, 22.05203933394581, 22.072336253551562, 22.092633173157314, 22.11293009276307, 22.13322701236882, 22.153523931974572, 22.173820851580324, 22.19411777118608, 22.254408518745116, 22.294402346698398, 22.334396174651683, 22.374390002604965, 22.41438383055825, 22.454377658511532, 22.494371486464818, 22.5343653144181, 22.574359142371385, 22.614352970324667, 22.67405524802451, 22.693763697771065, 22.713472147517624, 22.73318059726418, 22.752889047010736, 22.772597496757292, 22.79230594650385, 22.812014396250405, 22.831722845996964, 22.85143129574352, 22.967004935971776, 23.062870126453475, 23.158735316935175, 23.254600507416875, 23.350465697898574, 23.446330888380277, 23.542196078861977, 23.638061269343677, 23.733926459825376, 23.829791650307076, 23.981199219310742, 24.03674159783271, 24.092283976354675, 24.147826354876642, 24.20336873339861, 24.258911111920572, 24.31445349044254, 24.369995868964505, 24.42553824748647, 24.48108062600844, 24.578748185936732, 24.62087336734306, 24.662998548749385, 24.705123730155712, 24.74724891156204, 24.789374092968366, 24.831499274374693, 24.87362445578102, 24.915749637187346, 24.957874818593673, 25.0]
    
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

                    wd = 0.1
                    fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'
                    g = zeros(length(simul.y), length(simul.x))       

                    
                    

                    println("$omega: $Kx by $Ky")
                    
                    sol1, sol2, nIter = SEM_Wave_2d_Updated.Waveholtz(simul, omega, fVals, alphas, g, systemTol)
                    
                    

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
    
    println("values of omega: ")
    println(omegas)
    println("iterations for convergence: whi")
    println(whIters)
    println("iterations for convergence: gmres")
    println(gmresIters)


        
    # plot the result

    plt = plot(omegas, whIters, yscale=:log10, label="Waveholtz iteration", legend=:bottomleft)
    plot!(omegas, gmresIters, yscale=:log10, label="gmres iterations", legend=:bottomleft)

    savefig(plt, "BoxIterPlot")



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
    xr = 1.0
    yl = 0.0
    yr = 1.0

    alphas = [0.0; 0.0; 0.0; 0.0]
    
    omega = 5.0

    #spaceTol = 1e-3
    systemTol = 1e-12

    nIter = 5000

    N = 6   # degree of interpolation polynomials

    c_square(x, y) = 1 + 0.2*sin(x).*cos(x)

    #=
    L = maximum([xr-xl, yr-yl]) # length of the domain

    lambda = 1/(omega/(2*pi))
    N_lambda = L / lambda
    
    # rule of thumb 
    PPW = pi*(N_lambda/spaceTol)^(1/N)

    dofs_x = ceil(PPW * N_lambda)       # how many points do we want per dimension? points = point/wavelength * wavelength


    Kx = Int(ceil((dofs_x - 1)/N)) # the corresponding number of elements
    Ky = Kx # the same number in both directions!
    =#

    Kx = 4
    Ky = 4

    #Kx = 8
    #Ky = 8

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
    wd = 0.2
    fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'
    g = zeros(length(simul.y), length(simul.x))       

    ##################### byt till direkt lösare!!! #####################


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
    diffs_0_L2 = zeros(nIter, 1)
    diffs_1_L2 = zeros(nIter, 1)

    next_vs_now_0_H1 = zeros(nIter, 1)
    next_vs_now_1_H1 = zeros(nIter, 1)
    next_vs_now_seminorm = zeros(nIter, 1)

    # if we want inexact timestepping!
    simul.exactStep = false

    for j = 1:nIter

        uStart
        # perform WHI and measure H_1-differences with reference
        SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

        diffs_0[j] = SEM_Wave_2d_Updated.SobolevNorm(simul, simul.uFiltered - ref_0, 1.0)
        diffs_1[j] = SEM_Wave_2d_Updated.SobolevNorm(simul, simul.uDerFiltered - omega * ref_1, 1.0)
        diffs[j] = SEM_Wave_2d_Updated.SeminormWH(simul, simul.uFiltered - ref_0, simul.uDerFiltered - ref_1)

        # we also want to know how the quantities v_0^{n+1}-v_0^n and similar behave, so we measure these as well.
        next_vs_now_0_H1[j] = SEM_Wave_2d_Updated.SobolevNorm(simul, simul.uFiltered - uStart, 1.0)
        next_vs_now_1_H1[j] = SEM_Wave_2d_Updated.SobolevNorm(simul, simul.uDerFiltered - uStartDer, 1.0)
        next_vs_now_seminorm[j] = SEM_Wave_2d_Updated.SeminormWH(simul, simul.uFiltered - uStart, simul.uDerFiltered - uStartDer)


        # L2 norms may also be of interest
        diffs_0_L2[j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uFiltered, ref_0, 2)
        diffs_1_L2[j] = SEM_Wave_2d_Updated.LpNorm(simul, simul.uDerFiltered, omega * ref_1, 2)


        uStart .= simul.uFiltered
        uStartDer .= simul.uDerFiltered 

        println("iter $j of $nIter")
        
    end
                        
    # plot the result

    plt = plot(LinRange(1, nIter, nIter), diffs_0.^2 + diffs_1.^2, yscale=:log10, label="total iteration error", legend=:bottomleft)
    plot!(LinRange(1, nIter, nIter), diffs_0, yscale=:log10, label="v_0 error in H1")
    plot!(LinRange(1, nIter, nIter), diffs_1, yscale=:log10, label="v_1 error in H1")
    plot!(LinRange(1, nIter, nIter), diffs, yscale=:log10, label="seminorm error")



    savefig(plt, "BoxConvHistory")

    plt = heatmap(simul.x, simul.y, fVals)
    savefig(plt, "BoxRealPart.pdf")

    #=
    println("v_0 error in H1: ")
    println(diffs_0)
    println("v_1 error in H1: ")
    println(diffs_1)
    println("seminorm error: ")
    println(diffs)
    =#

    matwrite("BoxConvHistoryNewLong.mat", Dict(
    "v_0_H1" => diffs_0,
    "v_1_H1" => diffs_1,
    "v_0_L2" => diffs_0_L2,
    "v_1_L2" => diffs_1_L2,
    "seminorms" => diffs,
    "seq_v_0_H1" => next_vs_now_0_H1,
    "seq_v_1_H1" => next_vs_now_1_H1,
    "seq_seminorms" => next_vs_now_seminorm,
    "Re_U" => ref_0,
    "Im_U" => omega*ref_1,
    "v_0" => uStart,
    "v_1" => uStartDer,
    "xVals" => simul.x,
    "yVals" => simul.y,
    "problemData" => [omega, Kx, Ky, N] 
    ))


end




function convHistoryFinal()

    ###### setup ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    alphas = [0.0; 0.0; 0.0; 0.0]
    
    
    # setup for boxConV
    omega = 5.0
    nIter = 50
    N = 6   # degree of interpolation polynomials
    c_square(x, y) = 1
    

    
    # the data we want: second column has the data from the finer spatial discretisation
    e_0_L2 = zeros(nIter, 2)
    e_1_L2 = zeros(nIter, 2)
    
    grad_e_0_L2 = zeros(nIter, 2)
    grad_e_1_L2 = zeros(nIter, 2)

    # the same thing but with inexact timestepping
    e_0_L2_inexact = zeros(nIter, 2)
    e_1_L2_inexact = zeros(nIter, 2)
    grad_e_0_L2_inexact = zeros(nIter, 2)
    grad_e_1_L2_inexact = zeros(nIter, 2)



    number_of_refinements = 2


    for j = 1:number_of_refinements


        # choose number of elements
        Kx = Integer(4*j)
        Ky = Integer(4*j)


        # create simul object and forcing
        simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
        wd = 0.2
        
        fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'


        g = zeros(length(simul.y), length(simul.x))       

        # reference solution
        H, F = SEM_Wave_2d_Updated.HelmholtzMatrix(simul, omega, alphas, fVals)
        uSolVec = H\F
        uSol = reshape(uSolVec, length(simul.y), length(simul.x))
        u_re = real(uSol)
        u_im = imag(uSol)

        


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

            # note that GradIntegral computes the integral of c^2 \nabla U \cdot \nabla V
            grad_e_0_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uFiltered - u_re, simul.uFiltered - u_re))^(0.5)
            grad_e_1_L2[i, j] = (SEM_Wave_2d_Updated.GradIntegral(simul, simul.uDerFiltered - omega*u_im, simul.uDerFiltered - omega*u_im))^(0.5)

            uStart .= simul.uFiltered
            uStartDer .= simul.uDerFiltered 

            println("iter $i of $nIter, $Kx by $Kx elements, exact timestepping")
            
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

            uStart .= simul.uFiltered
            uStartDer .= simul.uDerFiltered 

            println("iter $i of $nIter, $Kx by $Kx elements, inexact timestepping")
            
        end

    end    



    matwrite("BoxConvHistoryAgain.mat", Dict(
    "e_0_L2" => e_0_L2,
    "e_1_L2" => e_1_L2,
    "grad_e_0_L2" => grad_e_0_L2,
    "grad_e_1_L2" => grad_e_1_L2,
    "e_0_L2_inexact" => e_0_L2_inexact,
    "e_1_L2_inexact" => e_1_L2_inexact,
    "grad_e_0_L2_inexact" => grad_e_0_L2_inexact,
    "grad_e_1_L2_inexact" => grad_e_1_L2_inexact
    ))


end

function convHistoryBoxNewNorm()


    # box model
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0


    alphas = [0.0; 0.0; 0.0; 0.0]
    
    omega = 5.0

    nIter = 1000

    N = 6   # degree of interpolation polynomials

    c_square(x, y) = 1


    KxBase = 4
    KyBase = 4
    
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

        
        wd = 0.2
        fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'

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



    matwrite("BoxConvHistoryNewNorm.mat", Dict(
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



function convHistoryVar()

    ###### setup ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    alphas = [0.0; 0.0; 0.0; 0.0]
    
    # a test with variable wavespeed
    omega = 5.0
    nIter = 25
    N = 6   # degree of interpolation polynomials
    # make bump function:
    g1(x) = exp(-1 ./ maximum([0, x]))
    a = 0.25; b = 0.35;
    #a = 0.499; b = 0.501;
    g2(x) = g1(x) ./ (g1(x) + g1(1-x))
    bump(x) = g2((x-a)/(b-a))
    bump2(x) = g2((x-0.75)/(0.85-0.75))

    #c_square(x, y) = 1.0 - 0.9*(bump(x) - bump2(x))
    c_square(x, y) = 1000*(1.0 - 0.5*bump(x))


    
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



    KxBase = 4
    KyBase = 4


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

    number_of_refinements = 2


    for j = 1:number_of_refinements


        # choose number of elements
        Kx = Integer(KxBase*j)
        Ky = Integer(KyBase*j)


        # create simul object and forcing
        simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)
        
        wd = 0.2
        fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'
        
        

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



    matwrite("BoxConvHistoryVar.mat", Dict(
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


#= the eigenvalues: 
   0.000000000000000
   3.141592653589793
   4.442882938158366
   6.283185307179586
   7.024814731040727
   8.885765876316732
   9.424777960769379
   9.934588265796101
  11.327173399138976
  12.566370614359172
  12.953118343415190
  13.328648814475097
  14.049629462081453
  15.707963267948966
  16.019042244414091
  16.917994196464051
  17.771531752633464
  18.318475636281480
  18.849555921538759
  19.109562078716149
  19.869176531592203
  20.116008064341784
  21.074444193122179
  21.991148575128552
  22.214414690791831
  22.654346798277953
  22.871139745490076
  23.925656840788776
  24.536623004530405
  25.132741228718345
  25.328329713402109
  25.906236686830379
  26.657297628950193
  26.841779398533230
  27.025001862730971
  28.099258924162907
  28.274333882308138
  28.448331425398703
  28.964053136475830
  29.637725818573745
  29.803764797388304
  30.941099316373158
  31.100180567108563
  31.415926535897931
  31.572615420804546
  32.038084488828183
  32.344676015002406
  32.799190229619086
  33.395587991875480
  33.835988392928101
  33.981520197416934
  35.124073655203631
  35.543063505266929
  35.819667392950713
  36.636951272562960
  37.829785066240554
  38.348025448024231
  39.985946443425291
  40.232016128683568
  42.265806470445746
  44.428829381583661

=#




function gmresSweepBox()
    
    ###### setup ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    
    alphas = [1/sqrt(2); 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)]
    #alphas = [0.0; 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)]
    #alphas = [0.0; 1/sqrt(2); 1/sqrt(2); 0.0]
    #alphas = [1/sqrt(2); 0.0; 1/sqrt(2); 0.0]
    #alphas = [1/sqrt(2); 0.0; 0.0; 0.0]
    #alphas = [0.0; 0.0; 0.0; 0.0]


    spaceTol = 1e-3
    systemTol = 1e-6

    
    h = 0.05
    omegaSmall = 2.0
    omegaBig = 15.0

    numberOfOmegas = Int(round((omegaBig-omegaSmall) / h + 1))
    omegas = collect(LinRange(omegaSmall, omegaBig, numberOfOmegas))

    #run for these values of omegas, record the number of iterations
    gmresIters = zeros(length(omegas))

    # iterate over the frequencies
    for j = 1:numberOfOmegas

        N = 6   # degree of interpolation polynomials

        c_square(x, y) = 1

        L = maximum([xr-xl, yr-yl]) # length of the domain

        lambda = 1/(omegas[j]/(2*pi))
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

        wd = 0.1
        fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'
        g = zeros(length(simul.y), length(simul.x))       
        

        seminormRes = Inf
        u_0_old = zeros(length(simul.y), length(simul.x))
        u_1_old = zeros(length(simul.y), length(simul.x))

        nIter = 0

        while (seminormRes > systemTol)

            nIter = nIter + 1            
            
            u_0, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omegas[j], fVals, alphas, g, nIter)

            seminormRes = SEM_Wave_2d_Updated.SeminormWH(simul, u_0 - u_0_old, u_1 - u_1_old) / SEM_Wave_2d_Updated.SeminormWH(simul, u_0_old, u_1_old)

            println("omega : $(omegas[j]) || nIter : $nIter || res : $seminormRes")

            u_0_old = u_0
            u_1_old = u_1

        end
        
        gmresIters[j] = nIter

    end

    println(omegas)
    println(gmresIters)

end



function makeBoxSolPlots()
    
    ###### setup ######
    
    xl = 0.0
    xr = 1.0
    yl = 0.0
    yr = 1.0

    omega = 12.75
    
    alphas_1 = [0.0; 1/sqrt(2); 1/sqrt(2); 1/sqrt(2)] # Neumann right
    alphas_2_adj = [0.0; 1/sqrt(2); 1/sqrt(2); 0.0] # Neumann right + down
    alphas_2_opp = [1/sqrt(2); 0.0; 1/sqrt(2); 0.0] # Neumann up + down
    alphas_3 = [1/sqrt(2); 0.0; 0.0; 0.0] # Neumann up + down + left
    alphas_4 = [0.0; 0.0; 0.0; 0.0] # Neumann all sides

    nIter_1 = 24
    nIter_2_adj = 30
    nIter_2_opp = 49
    nIter_3 = 73
    nIter_4 = 58

    Kx = 5; Ky = 5
    
    N = 6   # degree of interpolation polynomials
    c_square(x, y) = 1 # uniform wavespeed

    simul = SEM_Wave_2d_Updated.SEM_Wave([xl, xr], [yl, yr], Kx, Ky, N, c_square)

    # where we will store real parts of approximate solutions
    u_re_1 = zeros(length(simul.y), length(simul.x))
    u_re_2_adj = zeros(length(simul.y), length(simul.x))
    u_re_2_opp = zeros(length(simul.y), length(simul.x))
    u_re_3 = zeros(length(simul.y), length(simul.x))
    u_re_4 = zeros(length(simul.y), length(simul.x))

    wd = 0.1
    fVals = 10*exp.(-((simul.y .- 0.223)/wd).^2) * exp.(-((simul.x .- 0.652)/wd).^2)'
    g = zeros(length(simul.y), length(simul.x))       

    
    u_re_1, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas_1, g, nIter_1)
    u_re_2_adj, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas_2_adj, g, nIter_2_adj)
    u_re_2_opp, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas_2_opp, g, nIter_2_opp)
    u_re_3, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas_3, g, nIter_3)
    u_re_4, u_1, history = SEM_Wave_2d_Updated.WaveholtzGMRES(simul, omega, fVals, alphas_4, g, nIter_4)


    matwrite("BoxSolPlots.mat", Dict(
    "xVals" => simul.x,
    "yVals" => simul.y,
    "problemData" => [omega, Kx, Ky, N],
    "fVals" => fVals,
    "u_re_1" => u_re_1,
    "u_re_2_adj" => u_re_2_adj,
    "u_re_2_opp" => u_re_2_opp,
    "u_re_3" => u_re_3,
    "u_re_4" => u_re_4
    ))

end