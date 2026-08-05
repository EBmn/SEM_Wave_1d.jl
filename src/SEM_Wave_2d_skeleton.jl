function Simulate(simul::SEM_Wave, uStart, uStartDer, Tend::Float64, nsteps::Int64, forcing::Matrix{Float64}, omega::Float64, alphas::Vector{Float64}, g::Matrix{Float64}, animate = false, snapshotFrequency = 10, plotHeight = 1000.0, animationName = "test.gif")
    # approximates a solution to the wave equation using SEM_Wave
    # creates a gif named animationName if animate == true
    # snapshotFrequency decides how many steps are made before a frame is saved in the gif
    # the y-axis shown in the gif is [-plotHeight, plotHeight]

    # the filtered solution and derivative should always be reset, otherwise repeated calls to Waveholtz will produce the wrong result
    #to = TimerOutput()
    
    Initialise!(simul, uStart, uStartDer, Tend, nsteps, forcing, omega, alphas, g)  # calculates values for t = -simul.timestep using values at t = 0:


    # create an Animation if one is asked for
    if (animate == true)
       
        anim = Animation()

        surface(simul.x, simul.y[end:-1:1], simul.uNow, zlims = (-plotHeight, plotHeight), legend=:false)
        frame(anim)    

    end
    
    
    # make the steps

    for n = 1:simul.nsteps

        MakeStep!(simul, n)    

        # making the Animation is very slow compared to the actual calculations made here.
        

        # store a frame in our Animation once every snapshotFrequency steps
        if (animate == true && n%snapshotFrequency == 0)
            surface(simul.x, simul.y[end:-1:1], simul.uNow, zlims = (-plotHeight, plotHeight), legend=:false)
            frame(anim)
        end

            
    end

    # finally make and store the gif
    if (animate == true)

        gif(anim, animationName, fps=40) # if you need the gif to run a different framerate, change the value of fps

    end

    # return simul.uNow

    #show(to)

end



function MakeStep!(simul::SEM_Wave, stepnumber::Int64) 
    
    to = TimerOutput()

    t = stepnumber*simul.timestep # the time at uNow after(!) this step

    
    
    simul.uNext = TimeSteppingTerm(simul) + 
                  LaplaceTerm(simul, simul.uNow) + 
                  ForcingTerm(simul, stepnumber) + 
                  BoundaryTerm(simul, stepnumber)
    

                  
    # update the variables
    simul.uPrev = simul.uNow
    simul.uNow = simul.uNext

end




function LaplaceTerm(simul::SEM_Wave, U::Matrix{Float64})

    to = TimerOutput()

    @timeit to "fetching data" pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        

    @timeit to "make" laplaceVals = zeros(length(simul.y), length(simul.x))

    @timeit to "make" U_k = zeros(pointsPerElement, pointsPerElement)   
    @timeit to "make" V_k = zeros(pointsPerElement, pointsPerElement)

    @timeit to "make" V_kx = zeros(pointsPerElement, pointsPerElement) # x-derivative part
    @timeit to "make" V_ky = zeros(pointsPerElement, pointsPerElement) # y-derivative part
    @timeit to "make" scratch = zeros(pointsPerElement, 1)


    @timeit to "fetching data" G = simul.G
    @timeit to "fetching data" D = simul.D

    @timeit to "make" scratchMat = zeros(pointsPerElement, pointsPerElement)

    for k = 1:Kx*Ky
       
        @timeit to "get" U_k .= GetDegreesOfFreedom(simul, k, U)

        #@timeit to "get" c_square_k = GetDegreesOfFreedom(simul, k, simul.c_square) no longer needed!

        for j = 1:pointsPerElement # fill the matrices, column by column or row by row:

            @timeit to "VecMat" mul!((@view V_kx[j, :]), simul.DxMatrices[(k-1)*pointsPerElement + j], @view U_k[j, :])
            @timeit to "MatVec" mul!((@view V_ky[:, j]), simul.DyMatrices[(k-1)*pointsPerElement + j], @view U_k[:, j])

        end

    
        @timeit to "weightsMul" V_kx .= simul.QuadWeights .* V_kx
        @timeit to "weightsMul" V_ky .= simul.QuadWeights' .* V_ky
        
       
        @timeit to "add" V_k .= V_kx .+ V_ky
        
        @timeit to "set" SetDegreesOfFreedom!(simul, k, laplaceVals, V_k, true)
       
    end

    a = simul.stepCoeffs[1]; b = simul.stepCoeffs[2]

    @timeit to "adjustments" laplaceVals .= .- (simul.timestep^2 ./ ((simul.M ./ a) .+ (1 ./(2 .*b)) .*simul.timestep .* simul.bcMat .* simul.M_b)) .* laplaceVals   #u Laplace v = div(u grad v) - grad u \cdot grad v, hence the sign

    #show(to)

    return laplaceVals

end





function ForcingTerm(simul::SEM_Wave, stepnumber::Int64)

    pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        

    forcingVals = zeros(length(simul.y), length(simul.x))
    F_k = zeros(pointsPerElement, pointsPerElement)

    QuadWeightsMat = simul.QuadWeights * simul.QuadWeights'

    t = simul.timestep*(stepnumber-1)

    if simul.useMMS
        for i = 1:length(simul.y)
            for j = 1:length(simul.x)
                simul.fVals[i, j] = simul.c_square[i, j] * MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 2, 0, 0, simul.MMS_j) + # Laplacian!
                                    simul.c_square[i, j] * MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 0, 2, 0, simul.MMS_j) + 
                                    MMS.MMSfun(simul.x[j], simul.y[end+1-i], 1, 0, simul.MMS_c) * MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 1, 0, 0, simul.MMS_j) + 
                                    MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0, 1, simul.MMS_c) * MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 0, 1, 0, simul.MMS_j) - 
                                    MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 0, 0, 2, simul.MMS_j) # second-derivative in time
            end
        end
    end


    for k = 1:Kx*Ky

        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = simul.xNodes[i + 1] - simul.xNodes[i]
        delta_y_k = simul.yNodes[j + 1] - simul.yNodes[j]

        F_k .= 0.25*delta_x_k*delta_y_k * QuadWeightsMat .* GetDegreesOfFreedom(simul, k, simul.fVals) 

        SetDegreesOfFreedom!(simul, k, forcingVals, F_k, true)

    end

    if (simul.useMMS == false)
        forcingVals .= forcingVals .* cos(simul.omega*t)
    end

    #println("newer version:")
    #println(-(simul.timestep^2 ./ (simul.M + 0.5*simul.timestep * simul.bcMat .* simul.M_b)) .* forcingVals)
    a = simul.stepCoeffs[1]; b = simul.stepCoeffs[2]

    return -(simul.timestep^2 ./ ((simul.M ./ a) .+ (1 ./(2 .* b)) .* simul.timestep .* simul.bcMat .* simul.M_b)) .* forcingVals    # negative sign because we are solving u_tt = u_xx + u_yy - fcos(\omega x) as opposed to u_tt = u_xx + u_yy + fcos(\omega x).

end



function BoundaryTerm(simul::SEM_Wave, stepnumber::Int64)

    boundaryVals = zeros(length(simul.y), length(simul.x)) # not efficient, most entries of the matrix will be zero

    t = simul.timestep*(stepnumber-1) # current timepoint

    #alpha = simul.bc[1]; beta = simul.bc[2] # the case where beta = 0 has to be handled separately...
    QuadWeights = simul.QuadWeights
    xNodes = simul.xNodes
    yNodes = simul.yNodes

    pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        




    # perhaps store this rather than allocating again every time... Also: corners?
    betaInvMat = spzeros(length(simul.y), length(simul.x))
    
    # first attempt:
    
    
    betaInvMat[:, end] .+= 1/sqrt(1 - simul.alphas[1]^2)
    betaInvMat[1, :] .+= 1/sqrt(1 - simul.alphas[2]^2)
    betaInvMat[:, 1] .+= 1/sqrt(1 - simul.alphas[3]^2)
    betaInvMat[end, :] .+= 1/sqrt(1 - simul.alphas[4]^2)
    
    


    
    #=
    # second attempt: 
    # first add all betas...
    betaInvMat[:, end] .+= sqrt(1 - simul.alphas[1]^2)
    betaInvMat[1, :] .+= sqrt(1 - simul.alphas[2]^2)
    betaInvMat[:, 1] .+= sqrt(1 - simul.alphas[3]^2)
    betaInvMat[end, :] .+= sqrt(1 - simul.alphas[4]^2)

    # ...then flip them
    betaInvMat[:, 1] = 1 ./ betaInvMat[end, :]
    betaInvMat[:, end] = 1 ./ betaInvMat[:, end]
    betaInvMat[2:end-1, 1] = 1 ./ betaInvMat[2:end-1, 1] # careful to not flip the corners twice
    betaInvMat[2:end-1, end] = 1 ./ betaInvMat[2:end-1, end]
    =#

    if simul.useMMS

        # if we use MMS, g should be found using the target solution
        #weightedRefNormalDer, refTimeDer = SEM_Wave_2d_Updated.BoundaryDersMMS(simul, t)

        refNormalDer, refTimeDer = SEM_Wave_2d_Updated.BoundaryDersMMS(simul, t)

        alphaMat = spzeros(length(simul.y), length(simul.x))
        alphaMat[:, end] .+= simul.alphas[1]
        alphaMat[1, :] .+= simul.alphas[2]
        alphaMat[:, 1] .+= simul.alphas[3]
        alphaMat[end, :] .+= simul.alphas[4]

        # perhaps average the values at corners?
        #=
        alphaMat[1, 1] = alphaMat[1, 1]/2; alphaMat[1, end] = alphaMat[1, end]/2; alphaMat[end, 1] = alphaMat[end, 1]/2; alphaMat[end, end] = alphaMat[end, end]/2
        betaMat[1, 1] = betaMat[1, 1]/2; betaMat[1, end] = betaMat[1, end]/2; betaMat[end, 1] = betaMat[end, 1]/2; betaMat[end, end] = betaMat[end, end]/2
        # does not fix the issue
        =# 

            

        simul.g = simul.bcMat .* refTimeDer + simul.c_square .* refNormalDer
        #simul.g = alphaMat .* refTimeDer + weightedRefNormalDer

        boundaryVals = simul.M_b .* simul.g

    else
        # if no MMS, simul.g already stores the correct information about g.
        # Now integrate g along the boundary
        boundaryVals = simul.M_b .* (betaInvMat .* simul.g)


    end


    a = simul.stepCoeffs[1]; b = simul.stepCoeffs[2]

    return (simul.timestep^2 ./ ((simul.M ./ a) .+ (1 ./ (2 .* b)) .* simul.timestep.* simul.bcMat .* simul.M_b)) .* boundaryVals 

end



function TimeSteppingTerm(simul::SEM_Wave)

    a = simul.stepCoeffs[1]; b = simul.stepCoeffs[2]
    term = (2/a) .* simul.M .* simul.uNow + ((1 ./(2*b)) .* simul.timestep .* simul.bcMat .* simul.M_b - (simul.M ./ a)) .* simul.uPrev

    #println("newer version: ")
    #println(term ./ (simul.M + 0.5 * simul.timestep * simul.bcMat .* simul.M_b))

    return term ./ ((simul.M ./ a) + (1 ./(2 .* b)) .* simul.timestep .* simul.bcMat .* simul.M_b)

end



function Initialise!(simul::SEM_Wave, uStart::Matrix{Float64}, uStartDer::Matrix{Float64}, Tend::Float64, nsteps::Int64, forcing::Matrix{Float64}, omega::Float64, alphas::Vector{Float64}, g::Matrix{Float64})

    #time-related things
    simul.Tend = Tend
    simul.nsteps = nsteps
    simul.timestep = Tend/nsteps

    simul.omega = omega
    simul.fVals = forcing

    simul.alphas = alphas
    simul.g = g

    # reset bcMat and fill it with the correct values. Note the corner business!
    SetBCmat!(simul, alphas)

    simul.uNow = uStart

    simul.stepCoeffs[1] = (sin(omega*simul.timestep/2)/(omega*simul.timestep/2))^2 # the value to do with u_tt, denoted by a in the notes
    simul.stepCoeffs[2] = sin(omega*simul.timestep)/(omega*simul.timestep)         # the value to do with u_t, denoted by b in the notes

    a = simul.stepCoeffs[1]
    b = simul.stepCoeffs[2]

    # perhaps if bcMat and simul.M_b always come together one might want to combine the two. Especially if we want alpha to vary with x, y sometime later...

    #SU = (simul.M + 0.5 * simul.timestep * simul.bcMat .* simul.M_b) .* LaplaceTerm(simul, uStart) / simul.timestep^2 without exact timestepping
    SU = ((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* LaplaceTerm(simul, uStart) / simul.timestep^2

    G = ((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* BoundaryTerm(simul, 1) / simul.timestep^2
    F = ((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* ForcingTerm(simul, 1) / simul.timestep^2

    # old attempt at exact timestepping
    #simul.uPrev = uStart - simul.timestep*uStartDer + 0.5 * simul.timestep^2 * (a ./ simul.M).*(G - simul.bcMat .* simul.M_b .* uStartDer +
    #                                                                                        SU +
    #                                                                                        F)

    # new attempt at exact timestepping
    simul.uPrev = uStart - b * simul.timestep*uStartDer + 0.5 * simul.timestep^2 * (a ./ simul.M).*(G - simul.bcMat .* simul.M_b .* uStartDer +
                                                                                            SU +
                                                                                            F)


end


function SetBCmat!(simul::SEM_Wave, alphas)


    simul.bcMat = spzeros(length(simul.y), length(simul.x))
    simul.bcMat[:, end] .+= alphas[1]/(sqrt(1 - alphas[1]^2))
    simul.bcMat[1, :] .+= alphas[2]/(sqrt(1 - alphas[2]^2))
    simul.bcMat[:, 1] .+= alphas[3]/(sqrt(1 - alphas[3]^2))
    simul.bcMat[end, :] .+= alphas[4]/(sqrt(1 - alphas[4]^2))

    simul.bcMat[1, end] = (alphas[1] + alphas[2])/(sqrt(1 - alphas[1]^2) + sqrt(1 - alphas[2]^2))
    simul.bcMat[1, 1] = (alphas[2] + alphas[3])/(sqrt(1 - alphas[2]^2) + sqrt(1 - alphas[3]^2))
    simul.bcMat[end, 1] = (alphas[3] + alphas[4])/(sqrt(1 - alphas[3]^2) + sqrt(1 - alphas[4]^2))
    simul.bcMat[1, end] = (alphas[4] + alphas[1])/(sqrt(1 - alphas[4]^2) + sqrt(1 - alphas[1]^2))

end