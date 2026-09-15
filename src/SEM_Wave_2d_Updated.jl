

module SEM_Wave_2d_Updated

using LinearAlgebra
using SparseArrays
using Plots
using FastGaussQuadrature
using Polynomials
using TimerOutputs
using IterativeSolvers
using LinearMaps
using StaticArrays
#using LoopVectorization
using LinearAlgebra.BLAS
using Base.Threads

include("MMS.jl")
using .MMS


export SEM_Wave

mutable struct SEM_Wave
    # stores data for SEM approximation of the wave equation \frac{\partial^2 u}{\partial t^2} = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} - f(x)\cos(\omega t)

    xNodes::Vector{Float64}                   # the endpoints of the elements in the x-direction
    yNodes::Vector{Float64}                   # the endpoints of the elements in the y-direction
    N::Int64                                  # the degree of the Lagrange polynomials 
    x::Vector{Float64}                        # the actual coordinates between xNodes[1] and xNodes[end] where we approximate our solution
    y::Vector{Float64}                        # the actual coordinates between yNodes[1] and yNodes[end] where we approximate our solution
    c_square::Matrix{Float64}                 # matrix storing the value of the squared wave speed at each point in simul.x * simul.y'

    nsteps::Int64                             # number of steps we make in time
    Tend::Float64                             # the endpoint of our interval in time
    timestep::Float64                         # size of the timesteps (= Tend/nsteps)
    
    uPrev::Matrix{Float64}                    # value of the displacement in the previous timepoint
    uNow::Matrix{Float64}                     # value of the displacement in the current timepoint
    uNext::Matrix{Float64}                    # value of the displacement in the next timepoint
    scratch::Matrix{Float64}                  # scratch space
    uFiltered::Matrix{Float64}                # the filtered time average of u for use in Waveholtz
    uDerFiltered::Matrix{Float64}             # the filtered time average of u_t for use in Waveholtz

    fVals::Matrix{Float64}                    # values of the driving term of the equation
    omega::Float64                            # frequency parameter omega
    alphas::Vector{Float64}                   # value of alpha on the different boundaries: note that this needs to be constant along each boundary. Stored as [alpha_right; alpha_up; alpha_left; alpha_down]
    bcMat::SparseMatrixCSC{Float64, Int64}    # stores the alpha/beta-quotients that are used repeatedly in the method.
    g::Matrix{Float64}                        # function values for g in the boundary condition
    stepCoeffs::Vector{Float64}               # vector which stores the coefficients that are used for the exact timestepping 


    QuadPoints::Vector{Float64}               # the points in [-1, 1] we use in our Gauss-Lobatto quadrature
    QuadWeights::Vector{Float64}              # the weights corresponding to the points in QuadPoints
    D::Array{Float64}                         # matrix storing derivatives of Lagrange polynomials: D_{ij} = l'_j(\xi_i)
    G::Array{Float64}                         # the reference stiffness matrix: G_{jm} = \sum_{l=0}^{N} w_l l'_m(\xi_l) l'_j(\xi_l); w_l \in QuadWeights, \xi_l \in QuadPoints, l_j's are the Lagrange polynomials                            
    M::Matrix{Float64}                        # the mass matrix
    M_b::Matrix{Float64}                      # the mass matrix on the boundary, used to handle terms from the boundary intergral
    DxMatrices::Vector{Matrix{Float64}}
    DyMatrices::Vector{Matrix{Float64}}
    
    useMMS::Bool                              # determines whether an MMS test is run or not
    exactStep::Bool                           # determines whether exact timestepping is used
    MMS_j::MMS_jet                            # stores the relevant data for the MMS test
    MMS_c::MMS_jet                            # stores the relevant data for the MMS test with variable wave speed
    MMS_intDiff::Matrix{Float64}              # the integral of the absolute difference between the approximation and the MMS-reference at each point
    MMS_int::Matrix{Float64}                  # the integral of the absolute value of the solution in time, for use for computing relative errors


    function SEM_Wave(
        xLims::Vector{Float64},               # smallest and largest x-values
        yLims::Vector{Float64},               # smallest and largest y-values
        Kx::Int64,                            # number of elements in x-direction
        Ky::Int64,                            # number of elements in y-direction
        N::Int64,                             # degree of polynomials; every element has N + 1 points
        c_square::Function                    # wave speed
        )
        
        QuadPoints, QuadWeights = gausslobatto(N + 1)
        uPrev = zeros(Ky*N + 1, Kx*N + 1)
        uNow = zeros(Ky*N + 1, Kx*N + 1)
        uNext = zeros(Ky*N + 1, Kx*N + 1)
        scratch = zeros(Ky*N + 1, Kx*N + 1)
        
        xNodes = collect(LinRange(xLims[1], xLims[2], Kx+1))
        yNodes = collect(LinRange(yLims[1], yLims[2], Ky+1))
        x = ConstructX(xNodes, QuadPoints)
        y = ConstructX(yNodes, QuadPoints)
        c_squareMatrix = c_square.(x', y) # store values of the wave speed squared 


        fVals = zeros(Ky*N + 1, Kx*N + 1)
        omega = 0.0

        alphas = [0; 0; 0; 0] # Neumann by default.
        bcMat = spzeros(length(y), length(x)) # happens to be the correct values for Neumann bc

        g = zeros(Ky*N + 1, Kx*N + 1)

        uFiltered = zeros(Ky*N + 1, Kx*N + 1)
        uDerFiltered = zeros(Ky*N + 1, Kx*N + 1)
        

        useMMS = false
        exactStep = true        # we should use exact timestepping by default!

        MMS_intDiff = zeros(Ky*N + 1, Kx*N + 1)
        MMS_int = zeros(Ky*N + 1, Kx*N + 1)

        Tend = 1                # arbitrary values, just here as placeholders
        nsteps = 100
        timestep = Tend/nsteps
        stepCoeffs = ones(2)

        # make G and the inverse of the mass matrix
        G = ConstructG(QuadPoints, QuadWeights)
        D = ConstructD(QuadPoints, QuadWeights)
        M, M_b = ConstructMs(xNodes, yNodes, QuadWeights)
        M_b = sqrt.(c_squareMatrix) .* M_b

        # make and store the many small matrices that will be used to calculate the derivatives in LaplaceTerm.
        DxMatrices, DyMatrices = ConstructDmatrices(Kx, Ky, D, QuadWeights, c_squareMatrix, xNodes, yNodes) 

        # change some coefficients in the MMS
        MMS_j = MMS.MMS_jet(ones(3, 1), 3) # dimension 3, corresponding to (x, y, t) in that order.
        MMS_j.coeff[1, 2] = 2.0
        MMS_j.coeff[2, 1] = exp(1)
        MMS_j.coeff[2, 3] = -2
        MMS_j.coeff[1, 3] = (1 + sqrt(5))/2

        MMS_c = MMS.MMS_jet(ones(2, 1), 2) # dimension 2, corresponding to (x, y) in that order.

            new(xNodes, yNodes, N, x, y, c_squareMatrix, nsteps, Tend, timestep,
                uPrev, uNow, uNext, scratch, uFiltered, uDerFiltered, fVals, omega, alphas, bcMat, g, stepCoeffs, QuadPoints, QuadWeights, D, G, M, M_b, DxMatrices, DyMatrices, useMMS, exactStep, MMS_j, MMS_c, MMS_intDiff, MMS_int)
        
    end

end

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

        #println(norm(simul.uNow - uStart*cos(omega*n*(Tend/nsteps)))/norm(simul.uNow))

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


    
    #println("time: $t, cos(t):, $(cos(simul.omega*t)), value: $(simul.uNext[14, 15])")
    
    
    #println(maximum(abs.(LaplaceTermOld(simul, simul.uNow) - LaplaceTerm(simul, simul.uNow))))
    #println(maximum(abs.(ParLaplaceTerm(simul, simul.uNow) - LaplaceTerm(simul, simul.uNow))))

    #@timeit to "laplace" LaplaceTerm(simul, simul.uNow)
    #@timeit to "parallel laplace" ParLaplaceTerm(simul, simul.uNow)
    #@timeit to "laplace" LaplaceTerm(simul, simul.uNow)
    #@timeit to "laplace old" LaplaceTermOld(simul, simul.uNow)

    #=
    @timeit to "timestep" TimeSteppingTerm(simul)
    @timeit to "laplace" LaplaceTerm(simul, simul.uNow)
    @timeit to "forcing" ForcingTerm(simul, stepnumber)
    @timeit to "boundary" BoundaryTerm(simul, stepnumber)
    =#


    #=
    ########################### TRAPEZOIDAL RULE WITH uDer MODIFIED FOR EXACT TIMESTEPPING ###########################
    
    # use a second-order approximations of the derivative at t-simul.timestep using known values
    # note how using a central difference requires us to take an additional timestep for all times relevant to WaveHoltz to be considered
    uDer = (simul.uNext - simul.uPrev) ./ (2*simul.timestep*simul.stepCoeffs[2])

    
    diff = uDer ./ simul.uNow .+ simul.omega*tan(simul.omega*(t-simul.timestep))
    relDiff = diff ./ (simul.omega*tan(simul.omega*(t-simul.timestep)))

    #println("should be zero if imag(V) is zero: " * string(maximum(abs.((diff)))))
    #println("should be zero if imag(V) is zero (relDiff): " * string(norm((relDiff))))
    #println("should be zero if imag(V) is zero: " * string(maximum(abs.((relDiff)))))

    #=
    kvot = simul.uNext ./ simul.uNow
    kvot_target = cos(simul.omega*t)/cos(simul.omega*(t-simul.timestep))
    println("should also be zero if imag(V) is zero: " * string(norm((kvot .- kvot_target)./kvot_target)))
    println("should also be zero if imag(V) is zero (relative): " * string(norm((kvot .- kvot_target)./kvot_target)))
    =#

    if stepnumber == 1 || stepnumber == simul.nsteps
        weight = simul.timestep/2
    else
        weight = simul.timestep
    end


    ###################################################################

    #println("time: " * string(t) * " || maximum: " * string(maximum(simul.uNow)))    

    # value of K(t) at the time corresponding to simul.uNow
    Kval = (cos(simul.omega * (t - simul.timestep)) - 0.25) * (2/simul.Tend)


    simul.uFiltered = simul.uFiltered + weight .* simul.uNow .* Kval
    simul.uDerFiltered = simul.uDerFiltered + weight .* uDer .* Kval

    #println("time: " * string(t) * " || filtered maximum: " * string(maximum(simul.uFiltered)))
    #println("time: " * string(t) * " || filtered der maximum: " * string(maximum(simul.uDerFiltered)))    

    ###################################################################

    =#

    #=
    ########################### OLD VERSION ###########################
    
    # use a second-order approximations of the derivative at t and t-simul.timestep using known values
    uDer = ((3/2)*simul.uNext - 2*simul.uNow + 0.5*simul.uPrev) ./ simul.timestep 
    uDerPrev = (simul.uNext - simul.uPrev) ./ (2*simul.timestep)


    ###########
    # update the variables
    simul.uPrev = simul.uNow
    simul.uNow = simul.uNext
    #println("time: " * string(t) * " || maximum: " * string(maximum(simul.uNow)))    


    # filter the solution for use in Waveholtz:
    Know = (cos(simul.omega*t) - 0.25) * (2/simul.Tend)
    Kprev = (cos(simul.omega * (t - simul.timestep)) - 0.25) * (2/simul.Tend)

    simul.uFiltered = simul.uFiltered + (simul.uPrev .* Kprev + simul.uNow .* Know) .* (simul.timestep/2)
    simul.uDerFiltered = simul.uDerFiltered + (uDerPrev .* Kprev + uDer .* Know) .* (simul.timestep/2)
    # note how this makes the 1 2 ... 2 1 pattern from the trapezoid method happen

    #println("time: " * string(t) * " || filtered maximum: " * string(maximum(simul.uFiltered)))
    #println("time: " * string(t) * " || filtered der maximum: " * string(maximum(simul.uDerFiltered)))    


    ###################################################################
=#    


#=    
    ########################### NEW VERSION ###########################
    

    #@timeit to "time derivative" uDer = (simul.uNext - simul.uPrev) ./ (2*simul.timestep) # derivative at time t - simul.timestep (second order)
    uDer = (simul.uNext - simul.uPrev) ./ (2*simul.timestep*simul.stepCoeffs[2]) # derivative at time t - simul.timestep (second order)

    # filter the solution for use in Waveholtz
    
    # choose the value of a_0 in the filter; ideally one would not redo this every iteration...
    if simul.exactStep
        a_0 = 0.25 - 0.25*(tan(pi/simul.nsteps))^2
    else
        a_0 = 0.25
    end

    K = (cos(simul.omega*(t-simul.timestep)) - a_0) * (2/simul.Tend)

    if (stepnumber == 1)
        weight = 0.5 * simul.timestep
    else
        weight = 1.0 * simul.timestep
    end

    @timeit to "trapezoidal rule real part" simul.uFiltered = simul.uFiltered + simul.uNow * K * weight # simul.uNow corresponds to time t - simul.timestep
    @timeit to "trapezoidal rule imag part" simul.uDerFiltered = simul.uDerFiltered + uDer * K * weight


    if (stepnumber == simul.nsteps) # add the final term, corresponding to time t
        
        #uDer = (simul.uNext - simul.uNow)/simul.timestep #first-order approximation of the derivative at t = simul.Tend
        @timeit to "time derivative" uDer = ((3/2)*simul.uNext - 2*simul.uNow + 0.5*simul.uPrev) ./ simul.timestep # derivative at time t (second order)

        K = (cos(simul.omega*t) - 0.25) * (2/simul.Tend) 

        @timeit to "trapezoidal rule real part" simul.uFiltered = simul.uFiltered + simul.uNext * K * 0.5 * simul.timestep # simul.uNext corresponds to time t = simul.Tend
        @timeit to "trapezoidal rule imag part" simul.uDerFiltered = simul.uDerFiltered + uDer * K * 0.5 * simul.timestep 

    end
    ###################################################################
=#

    



    #################### the left-Riemann sum that Amit uses ####################

    # second-order formulation => we need to approximate u_t
    #uDer = (simul.uNext - simul.uPrev) ./ (2*simul.timestep) # derivative at time t - simul.timestep (second order)

    uDer = (simul.uNext - simul.uPrev) ./ (2*simul.timestep*simul.stepCoeffs[2])

    # filter the solution for use in Waveholtz
    
    # choose the value of a_0 in the filter; for maximum elegance one would not redo this every iteration.
    if simul.exactStep
        a_0 = 0.25 - 0.25*(tan(pi/simul.nsteps))^2
    else
        a_0 = 0.25
    end

    K = (cos(simul.omega*(t-simul.timestep)) - a_0) * (2/simul.nsteps)

    simul.uDerFiltered = simul.uDerFiltered + uDer * K    
    simul.uFiltered = simul.uFiltered + simul.uNow * K # simul.uNow corresponds to time t - simul.timestep
    # we are using a left Riemann sum, so we do not use the data at time = T. (i.e. one should really not take that timestep in the first place)
    
    
    
    ###########################################################



    
    ###################################################################

    # update the variables
    simul.uPrev = simul.uNow
    simul.uNow = simul.uNext

    #println("time: $t")
    #println("max abs of u: " * string(maximum(abs.(simul.uNow))))


    # integrates the MMS error, if we are doing an MMS run.
    if (simul.useMMS)

        reference = zeros(length(simul.y), length(simul.x))

        # construct initial data and reference solution
        for i = 1:length(simul.y)
            for j = 1:length(simul.x)
                reference[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 0, 0, 0, simul.MMS_j)
            end
        end

        simul.MMS_intDiff += abs.(simul.uNow - reference) * 2 * simul.timestep
        simul.MMS_int += abs.(simul.uNow) * 2 * simul.timestep

        if (stepnumber == simul.nsteps) # do not count the final step twice!

            simul.MMS_intDiff -=(abs.(simul.uNow - reference)) * simul.timestep
            simul.MMS_int -= abs.(simul.uNow) * simul.timestep

        end

    end

    #show(to)


    

end

function LaplaceTermOld(simul::SEM_Wave, U::Matrix{Float64}) # this function is still here for archeological reasons, but not in use anywhere.

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

    alpha = simul.bc[1]; beta = simul.bc[2]




    # placeholders for doing the matmuls using mul! (probably incredibly suboptimal v.r.t. elegance)
    @timeit to "make" DtW = zeros(pointsPerElement, pointsPerElement) 
    @timeit to "make" prod1 = zeros(pointsPerElement, pointsPerElement) 
    @timeit to "make" prod2 = zeros(pointsPerElement, pointsPerElement) 
    @timeit to "make" row = zeros(1, pointsPerElement)
    @timeit to "make" column = zeros(pointsPerElement, 1)

    @timeit to "make" scratch1 = zeros(pointsPerElement, pointsPerElement) 
    @timeit to "make" scratch2 = zeros(pointsPerElement, pointsPerElement)

    @timeit to "diagm" W = diagm(simul.QuadWeights)

    @timeit to "mul clean" mul!(DtW, transpose(simul.D), W) # a matrix which is needed for both V_kx and V_ky

    for k = 1:Kx*Ky

        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = simul.xNodes[i+1] - simul.xNodes[i]
        delta_y_k = simul.yNodes[j+1] - simul.yNodes[j]
        
        @timeit to "get" U_k .= GetDegreesOfFreedom(simul, k, U)

        @timeit to "get" c_square_k = GetDegreesOfFreedom(simul, k, simul.c_square)

        
        
        @timeit to "make" V_kx = zeros(pointsPerElement, pointsPerElement) # x-derivative part
        @timeit to "make" V_ky = zeros(pointsPerElement, pointsPerElement) # y-derivative part

        # placeholders for doing the matmuls using mul! (probably incredibly suboptimal v.r.t. elegance)
        @timeit to "make" prod1 = zeros(pointsPerElement, pointsPerElement) 
        @timeit to "make" prod2 = zeros(pointsPerElement, pointsPerElement) 
        @timeit to "make" row = zeros(1, pointsPerElement)
        @timeit to "make" column = zeros(pointsPerElement, 1)

        for j = 1:pointsPerElement # fill the matrices, column by column or row by row:

            ########### older version ###########
            # build V_kx and V_xy in steps, adding a factor at a time...
            @timeit to "mul diag 1" mul!(prod1, Diagonal(c_square_k[j, :]), simul.D)
            @timeit to "mul clean" mul!(prod2, DtW, prod1)
            @timeit to "mul" mul!(row, Matrix(U_k[j, :]'), prod2, delta_y_k/delta_x_k, 0)
            @timeit to "store" V_kx[j, :] = row

            @timeit to "mul diag 1" mul!(prod1, Diagonal(c_square_k[:, j]), simul.D)
            @timeit to "mul clean" mul!(prod2, DtW, prod1)
            @timeit to "mul" mul!(column, prod2, U_k[:, j], delta_x_k/delta_y_k, 0)
            @timeit to "store" V_ky[:, j] = column

        end

    
        @timeit to "weightsMul" V_kx .= simul.QuadWeights .* V_kx
        @timeit to "weightsMul" V_ky .= simul.QuadWeights' .* V_ky
        
        #V_kx = (delta_y_k/delta_x_k) * simul.QuadWeights .* V_kx
        #V_ky = (delta_x_k/delta_y_k) * simul.QuadWeights' .* V_ky

        @timeit to "add" V_k .= V_kx .+ V_ky
        
        @timeit to "set" SetDegreesOfFreedom!(simul, k, laplaceVals, V_k, true)
        #SetDegreesOfFreedom!(simul, k, laplaceVals, V_k_old, true)

    end

    @timeit to "adjustments" laplaceVals .= - (simul.timestep^2 ./ (simul.M .+ 0.5 .* (alpha/beta).*simul.timestep.*simul.M_b)) .* laplaceVals   #u Laplace v = div(u grad v) - grad u \cdot grad v, hence the sign

    #show(to)

    return laplaceVals

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

    @timeit to "adjustments" laplaceVals .= .- (simul.timestep^2 ./ ((simul.M ./ a) .+ (simul.timestep ./(2 .* b)) .* simul.bcMat .* simul.M_b)) .* laplaceVals   #u Laplace v = div(u grad v) - grad u \cdot grad v, hence the sign

    #show(to)

    return laplaceVals

end



function ParLaplaceTerm(simul::SEM_Wave, U::Matrix{Float64}) ############# Claude-generated from LaplaceTerm #############, defunct

    to = TimerOutput()

    @timeit to "fetching data" pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        


    @timeit to "fetching data" G = simul.G
    @timeit to "fetching data" D = simul.D

    alpha = simul.bc[1]; beta = simul.bc[2]


    @timeit to "make" laplaceVals = zeros(length(simul.y), length(simul.x))

    #lk = ReentrantLock()

    used_buffers = Channel{Matrix{Float64}}(Inf)


    nbuffers = 2 * nthreads()
    buffers = [(
        U_k          = zeros(pointsPerElement, pointsPerElement),
        V_kx         = zeros(pointsPerElement, pointsPerElement),
        V_ky         = zeros(pointsPerElement, pointsPerElement),
        laplaceLocal = zeros(length(simul.y), length(simul.x)),
    ) for _ in 1:nbuffers]

    @timeit to "parallel laplace loop" @threads for k in 1:Kx*Ky
        
        buf = buffers[threadid()]  # no allocation, just an index lookup
        (; U_k, V_kx, V_ky, laplaceLocal) = buf

        # previous version
        #= 
        is_new = !haskey(task_local_storage(), :laplace_bufs)

        buf = get!(() -> (
            U_k     = zeros(pointsPerElement, pointsPerElement),
            V_kx    = zeros(pointsPerElement, pointsPerElement),
            V_ky    = zeros(pointsPerElement, pointsPerElement),
            laplaceLocal = zeros(length(simul.y), length(simul.x)),
        ), task_local_storage(), :laplace_bufs)

        if is_new
            put!(used_buffers, buf.laplaceLocal)
        end 

        =#

        #=
        # new attempt
        if !haskey(task_local_storage(), :laplace_bufs)
            task_local_storage(:laplace_bufs, (
                U_k          = zeros(pointsPerElement, pointsPerElement),
                V_kx         = zeros(pointsPerElement, pointsPerElement),
                V_ky         = zeros(pointsPerElement, pointsPerElement),
                laplaceLocal = zeros(length(simul.y), length(simul.x)),
            ))
            put!(used_buffers, task_local_storage(:laplace_bufs).laplaceLocal)
        end
        buf = task_local_storage(:laplace_bufs)

        (; U_k, V_kx, V_ky, laplaceLocal) = buf
        =#

        U_k .= GetDegreesOfFreedom(simul, k, U)

        for j = 1:pointsPerElement # fill the matrices, column by column or row by row:

            mul!((@view V_kx[j, :]), simul.DxMatrices[(k-1)*pointsPerElement + j], @view U_k[j, :])
            mul!((@view V_ky[:, j]), simul.DyMatrices[(k-1)*pointsPerElement + j], @view U_k[:, j])

        end


        V_kx .= simul.QuadWeights .* V_kx
        V_ky .= simul.QuadWeights' .* V_ky
        V_kx .+= V_ky

        SetDegreesOfFreedom!(simul, k, laplaceLocal, V_kx, true)
        

        #=
        lock(lk) do        
            SetDegreesOfFreedom!(simul, k, laplaceVals, V_kx, true)
        end
        =#

    end

    close(used_buffers)

    fill!(laplaceVals, 0)

    for buf in buffers
        laplaceVals .+= buf.laplaceLocal
    end



    #=
    @timeit to "merge thread data" for val in values(task_local_storage())
        if val isa NamedTuple && hasfield(typeof(val), :laplaceLocal)
            laplaceVals .+= val.laplaceLocal
        end
    end
    =#

    @timeit to "adjustments" laplaceVals .= .- (simul.timestep^2 ./ (simul.M .+ 0.5 .* (alpha/beta).*simul.timestep.*simul.M_b)) .* laplaceVals   #u Laplace v = div(u grad v) - grad u \cdot grad v, hence the sign

    show(to)

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

    #return -(simul.timestep^2 ./ ((simul.M ./ a) .+ (simul.timestep ./(2*b)) .* simul.bcMat .* simul.M_b)) .* (1/a) .* forcingVals    # negative sign because we are solving u_tt = u_xx + u_yy - fcos(\omega x) as opposed to u_tt = u_xx + u_yy + fcos(\omega x).
    return -(simul.timestep^2 ./ ((simul.M ./ a) .+ (simul.timestep ./(2*b)) .* simul.bcMat .* simul.M_b)) .* forcingVals    # negative sign because we are solving u_tt = u_xx + u_yy - fcos(\omega x) as opposed to u_tt = u_xx + u_yy + fcos(\omega x).

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




    

    

    #=
    term = (simul.timestep^2 ./ (simul.M + 0.5 * simul.timestep * simul.bcMat .* simul.M_b)) .* boundaryVals
    term[:, end] = term[:, end] ./ betaMat[:, end]
    term[1, 2:end-1] = term[1, 2:end-1] ./ betaMat[1, 2:end-1]
    term[:, 1] = term[:, 1] ./ betaMat[:, 1]
    term[end, 2:end-1] = term[end, 2:end-1] ./ betaMat[end, 2:end-1]
    =#

    #println("newer version!")
    #println((simul.timestep^2 ./ (simul.M + 0.5 * simul.timestep * simul.bcMat .* simul.M_b)) .* boundaryVals)    
    #println((simul.timestep^2 ./ (simul.M + 0.5 * simul.timestep * simul.bcMat .* simul.M_b)) .* boundaryVals)
    
    a = simul.stepCoeffs[1]; b = simul.stepCoeffs[2]

    return (simul.timestep^2 ./ ((simul.M ./ a) .+ (1 ./ (2 .* b)) .* simul.timestep.* simul.bcMat .* simul.M_b)) .* (1/b) .* boundaryVals 
    #return (1/beta) * (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* boundaryVals 

end


function TimeSteppingTerm(simul::SEM_Wave)

    a = simul.stepCoeffs[1]; b = simul.stepCoeffs[2]
    term = (2/a) .* simul.M .* simul.uNow + ((1 ./(2*b)) .* simul.timestep .* simul.bcMat .* simul.M_b - (simul.M ./ a)) .* simul.uPrev

    #println("newer version: ")
    #println(term ./ (simul.M + 0.5 * simul.timestep * simul.bcMat .* simul.M_b))

    return term ./ ((simul.M ./ a) + (1 ./(2 .* b)) .* simul.timestep .* simul.bcMat .* simul.M_b)

end


function TimeSteppingTerm!(simul::SEM_Wave) # OUT OF DATE: DO NOT USE YET!

    println("You are using an old TimeSteppingTerm function! Your output is probably incorrect!")

    simul.scratch .= 2*simul.M .* simul.uNow .+ (0.5 * simul.bcMat .* simul.timestep*simul.M_b - simul.M) .* simul.uPrev

    simul.uNext .= simul.uNext .+ simul.scratch ./ (simul.M + 0.5 * simul.bcMat * simul.timestep * simul.M_b)

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

    #=
    simul.bcMat = spzeros(length(simul.y), length(simul.x))
    simul.bcMat[:, end] .+= alphas[1]/(sqrt(1 - alphas[1]^2))
    simul.bcMat[1, :] .+= alphas[2]/(sqrt(1 - alphas[2]^2))
    simul.bcMat[:, 1] .+= alphas[3]/(sqrt(1 - alphas[3]^2))
    simul.bcMat[end, :] .+= alphas[4]/(sqrt(1 - alphas[4]^2))

    #simul.bcMat[1, end] = 0.5*simul.bcMat[1, end]
    #simul.bcMat[1, 1] = 0.5*simul.bcMat[1, 1]
    #simul.bcMat[end, 1] = 0.5*simul.bcMat[end, 1]
    #simul.bcMat[end, end] = 0.5*simul.bcMat[end, end]

    simul.bcMat[1, end] = (alphas[1] + alphas[2])/(sqrt(1 - alphas[1]^2) + sqrt(1 - alphas[2]^2))
    simul.bcMat[1, 1] = (alphas[2] + alphas[3])/(sqrt(1 - alphas[2]^2) + sqrt(1 - alphas[3]^2))
    simul.bcMat[end, 1] = (alphas[3] + alphas[4])/(sqrt(1 - alphas[3]^2) + sqrt(1 - alphas[4]^2))
    simul.bcMat[1, end] = (alphas[4] + alphas[1])/(sqrt(1 - alphas[4]^2) + sqrt(1 - alphas[1]^2))
    =#

    simul.uNow = uStart

    simul.uFiltered = zeros(length(simul.y), length(simul.x))
    simul.uDerFiltered = zeros(length(simul.y), length(simul.x))

    # if we are using MMS, we need to store the correct values in simul.c_square 

    if (simul.useMMS == true)

        for j = 1:length(simul.x)
            for i = 1:length(simul.y)

                simul.c_square[i, j] = 1.0 + MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0, 0, simul.MMS_c)

                # check that each value of c_square is reasonable
                if (simul.c_square[i, j] < 1e-6)
                    println("OBS: c_square is less than 1e-6 at point (" * string(simul.x[j]) * ", " * string(simul.y[end+1-i]) * ")")
                end    

            end

        end

        # we also need to update the precomputed matrices in that case
        DxMatrices, DyMatrices = ConstructDmatrices(length(simul.xNodes)-1, length(simul.yNodes)-1, simul.D, simul.QuadWeights, simul.c_square, simul.xNodes, simul.yNodes) 
        simul.DxMatrices = DxMatrices
        simul.DyMatrices = DyMatrices

    end

    if (simul.exactStep == true) # will happen unless the user very actively sets simul.exactStep = false

        simul.stepCoeffs[1] = (sin(0.5*omega*simul.timestep)/(0.5*omega*simul.timestep))^2 # the value to do with u_tt, denoted by a in the notes
        simul.stepCoeffs[2] = sin(omega*simul.timestep)/(omega*simul.timestep)         # the value to do with u_t, denoted by b in the notes

    else # if "normal" leapfrog timestepping should be used, that corresponds to values of 1.0. One could in principle put whichever value to get many kinds of schemes...

        simul.stepCoeffs[1] = 1.0
        simul.stepCoeffs[2] = 1.0

    end


    a = simul.stepCoeffs[1]
    b = simul.stepCoeffs[2]

    # perhaps if bcMat and simul.M_b always come together one might want to combine the two. Especially if we want alpha to vary with x, y sometime later...

    #SU = (simul.M + 0.5 * simul.timestep * simul.bcMat .* simul.M_b) .* LaplaceTerm(simul, uStart) / simul.timestep^2 without exact timestepping
    SU = ((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* LaplaceTerm(simul, uStart) / simul.timestep^2

    G = ((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* BoundaryTerm(simul, 1) / simul.timestep^2
    F = ((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* ForcingTerm(simul, 1) / simul.timestep^2


    # new attempt at exact timestepping
    simul.uPrev = uStart - b * simul.timestep*uStartDer + 0.5 * a * simul.timestep^2 * (1 ./ simul.M).*(G - simul.bcMat .* simul.M_b .* uStartDer +
                                                                                            SU + 
                                                                                            F)

    #simul.uPrev = uStart - b * simul.timestep*uStartDer + 0.5 * a * simul.timestep^2 * (1 ./ simul.M).*(F - simul.bcMat .* simul.M_b .* uStartDer + SU)

    #println("should be zero: " *string((1 - 0.5*simul.timestep^2 *omega^2*a) - cos(omega*simul.timestep)))
    #println("should be zero: " *string(norm((1 ./ simul.M).*(F - simul.bcMat .* simul.M_b .* uStartDer + SU) + simul.omega^2 * uStart)/norm(simul.omega^2 * uStart)))

    #=
    plt = heatmap(simul.y, simul.x, log10.(abs.((1 ./ simul.M).*(F - simul.bcMat .* simul.M_b .* uStartDer + SU) + simul.omega^2*uStart)))
    savefig(plt, "tester.pdf")
    plt = heatmap(simul.y, simul.x, (1 ./ simul.M).*(F - simul.bcMat .* simul.M_b .* uStartDer + SU))
    savefig(plt, "tester1.pdf")
    plt = heatmap(simul.y, simul.x, -simul.omega^2 * uStart)
    savefig(plt, "tester2.pdf")
    =#

    uPrevExpected = uStart*cos(simul.omega*simul.timestep) - (1/simul.omega) * uStartDer*sin(simul.omega*simul.timestep)
    #println("should be zero: " * string(norm(simul.uPrev - uPrevExpected)))

    if (simul.useMMS == true)
        for j = 1:length(simul.x)
            for i = 1:length(simul.y)
                simul.uPrev[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], -simul.timestep, 0, 0, 0, simul.MMS_j)
            end
        end
    end

end


function SetBCmat!(simul::SEM_Wave, alphas)


    simul.bcMat = spzeros(length(simul.y), length(simul.x))
    simul.bcMat[:, end] .+= alphas[1]/(sqrt(1 - alphas[1]^2))
    simul.bcMat[1, :] .+= alphas[2]/(sqrt(1 - alphas[2]^2))
    simul.bcMat[:, 1] .+= alphas[3]/(sqrt(1 - alphas[3]^2))
    simul.bcMat[end, :] .+= alphas[4]/(sqrt(1 - alphas[4]^2))

    #simul.bcMat[1, end] = 0.5*simul.bcMat[1, end]
    #simul.bcMat[1, 1] = 0.5*simul.bcMat[1, 1]
    #simul.bcMat[end, 1] = 0.5*simul.bcMat[end, 1]
    #simul.bcMat[end, end] = 0.5*simul.bcMat[end, end]

    simul.bcMat[1, end] = (alphas[1] + alphas[2])/(sqrt(1 - alphas[1]^2) + sqrt(1 - alphas[2]^2))
    simul.bcMat[1, 1] = (alphas[2] + alphas[3])/(sqrt(1 - alphas[2]^2) + sqrt(1 - alphas[3]^2))
    simul.bcMat[end, 1] = (alphas[3] + alphas[4])/(sqrt(1 - alphas[3]^2) + sqrt(1 - alphas[4]^2))
    simul.bcMat[1, end] = (alphas[4] + alphas[1])/(sqrt(1 - alphas[4]^2) + sqrt(1 - alphas[1]^2))

end


function ConstructX(nodes::Vector{Float64}, QuadPoints::Vector{Float64})
    
    K = length(nodes)-1
    N = length(QuadPoints)

    x = zeros(K*(N-1) + 1)

    # puts the correct values in x
    for k = 1:K
  
        delta_x_k = nodes[k+1]-nodes[k]

        for i = 1:N-1
            x[(k-1)*(N-1) + i] = nodes[k] + 0.5*(1 + QuadPoints[i])*delta_x_k 
        end

    end
    
    x[end] = nodes[end]

    return x

end


function ConstructG(QuadPoints::Vector{Float64}, QuadWeights::Vector{Float64})

    N = length(QuadWeights)

    D = zeros(N, N)

    weights = BarycentricWeights(QuadPoints)

    for i = 1:N
        for j = 1:N
            if (i != j)
                D[i, j] = (weights[j]/weights[i])*(1/(QuadPoints[i] - QuadPoints[j]))
                D[i, i] -= D[i, j]
            end
        end
    end

    return transpose(D)*(QuadWeights.*D)

end

function ConstructD(QuadPoints::Vector{Float64}, QuadWeights::Vector{Float64})

    N = length(QuadWeights)

    D = zeros(N, N)

    weights = BarycentricWeights(QuadPoints)

    for i = 1:N
        for j = 1:N
            if (i != j)
                D[i, j] = (weights[j]/weights[i])*(1/(QuadPoints[i] - QuadPoints[j]))
                D[i, i] -= D[i, j]
            end
        end
    end

    return D

end

function ConstructDmatrices(Kx, Ky, D, QuadWeights, c_squareMatrix, xNodes, yNodes)

    # we want to, for each element, create and store the matrices (delta_x_k/delta_y_k)D^T diag(QuadWeigths) diag(c^k[row j]) D
    # and the matrices (delta_y_k/delta_x_k) D^T diag(QuadWeights) diag(c^k[col j]) D
    # these are used in LaplaceTerm

    pointsPerElement = length(QuadWeights)
    #DxMatrices = Vector{SMatrix{pointsPerElement, pointsPerElement}}(undef, Kx*Ky*pointsPerElement)
    #DyMatrices = Vector{SMatrix{pointsPerElement, pointsPerElement}}(undef, Kx*Ky*pointsPerElement)
    DxMatrices = Vector{Matrix{Float64}}(undef, Kx*Ky*pointsPerElement)
    DyMatrices = Vector{Matrix{Float64}}(undef, Kx*Ky*pointsPerElement)
    W = spdiagm(QuadWeights)

    DtW = transpose(D)*W

    for k = 1:Kx*Ky
        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = xNodes[i+1] - xNodes[i]
        delta_y_k = yNodes[j+1] - yNodes[j]

        # practically a duplicate of GetDegreesOfFreedom. We cannot use GetDegreesOfFreedom as is, as we do not yet have a SEM_Wave object to pass into the function...
        xStart = (i - 1) * (pointsPerElement - 1)
        yStart = (j - 1) * (pointsPerElement - 1)

        c_square_k = @view c_squareMatrix[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)]
        

        for t = 1:pointsPerElement
            #DxMatrices[(k-1)*pointsPerElement + t] = SMatrix{pointsPerElement, pointsPerElement}((delta_y_k/delta_x_k) * DtW * spdiagm(c_square_k[t, :]) * D)
            #DyMatrices[(k-1)*pointsPerElement + t] = SMatrix{pointsPerElement, pointsPerElement}((delta_x_k/delta_y_k) * DtW * spdiagm(c_square_k[:, t]) * D)
            DxMatrices[(k-1)*pointsPerElement + t] = Matrix(transpose((delta_y_k/delta_x_k) * DtW * spdiagm(c_square_k[t, :]) * D)) # in fact we will really make use the transpose of that matrix
            DyMatrices[(k-1)*pointsPerElement + t] = (delta_x_k/delta_y_k) * DtW * spdiagm(c_square_k[:, t]) * D
        end

    end

    return DxMatrices, DyMatrices

end


function BarycentricWeights(points::Vector{Float64})

    # returns the barycentric weights of a vector of points; used in calculating the derivative matrix of the Lagrange polynomials

    N = length(points)
    weights = ones(N)

    for i = 1:N
        for j = 1:N
            if (i != j)
                weights[j] *= 1/(points[j] - points[i])
            end
        end
    end

    return weights

end


function ConstructMs(xNodes::Vector{Float64}, yNodes::Vector{Float64}, QuadWeights::Vector{Float64})

    # constructs the mass matrix and the boundary mass matrix

    pointsPerElement = length(QuadWeights)
    Kx = length(xNodes)-1
    Ky = length(yNodes)-1

    M = zeros(Ky*(pointsPerElement-1) + 1, Kx*(pointsPerElement-1) + 1)
    M_b = zeros(Ky*(pointsPerElement-1) + 1, Kx*(pointsPerElement-1) + 1)

    QuadWeightsMat = QuadWeights * QuadWeights'
    
    ### first construct the mass matrix M ###
    for k = 1 : Kx*Ky

        # "coordinates" of the element
        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1
            
        delta_x_k = xNodes[i+1] - xNodes[i]
        delta_y_k = yNodes[j+1] - yNodes[j]

        

        M_k = 0.25*delta_x_k * delta_y_k * QuadWeightsMat

        # practically a duplicate of SetDegreesOfFreedom. We cannot use SetDegreesOfFreedom as is, as we do not yet have a SEM_Wave object to pass into the function...
        xStart = (i - 1) * (pointsPerElement - 1)
        yStart = (j - 1) * (pointsPerElement - 1)

        M[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)] .+= M_k
        
    end

    ### next construct the boundary mass matrix M_b ###

    # M_b is mostly zeros and has nice symmetries

    # we go along the edges of the matrix
    
    for j = 1:Kx
        
        delta_x_k = xNodes[j+1] - xNodes[j]
        xStart = (j - 1) * (pointsPerElement - 1)

        M_b[end, (xStart + 1):(xStart + pointsPerElement)] .+= QuadWeights * delta_x_k * 0.5  # lower side
        M_b[1, (xStart + 1):(xStart + pointsPerElement)] .+= QuadWeights * delta_x_k * 0.5    # upper side

    end

    for j = 1:Ky

        delta_y_k = yNodes[end+1-j] - yNodes[end-j]
        yStart = (j - 1) * (pointsPerElement - 1)

        M_b[(yStart + 1):(yStart + pointsPerElement), 1] .+= QuadWeights * delta_y_k * 0.5    # left side
        M_b[(yStart + 1):(yStart + pointsPerElement), end] .+= QuadWeights * delta_y_k * 0.5  # right side

    end

    return M, M_b

end



function GetDegreesOfFreedom(simul::SEM_Wave, k::Int64, u::AbstractMatrix{<:Union{Float64,ComplexF64}})
    
    pointsPerElement = simul.N + 1
    Kx = length(simul.xNodes) - 1
    Ky = length(simul.yNodes) - 1

    # "coordinates" of the element
    i = mod(k - 1, Kx) + 1
    j = Int((k - i) / Kx) + 1

    xStart = (i - 1) * (pointsPerElement - 1)
    yStart = (j - 1) * (pointsPerElement - 1)

    return @view u[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)]
    #return u[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)]

end


function SetDegreesOfFreedom!(simul::SEM_Wave, k::Int64, v::AbstractMatrix{<:Union{Float64,ComplexF64}}, v_k::AbstractMatrix{<:Union{Float64,ComplexF64}}, add::Bool)

    pointsPerElement = simul.N + 1
    Kx = length(simul.xNodes) - 1
    Ky = length(simul.yNodes) - 1

    # "coordinates" of the element
    i = mod(k - 1, Kx) + 1
    j = Int((k - i) / Kx) + 1

    xStart = (i - 1) * (pointsPerElement - 1)
    yStart = (j - 1) * (pointsPerElement - 1)

    if add
        #v[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)] .+= v_k
        @views v[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)] .+= v_k
    else
        #v[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)] .= v_k
        @views v[(yStart + 1):(yStart + pointsPerElement), (xStart + 1):(xStart + pointsPerElement)] .= v_k
    end

end

########################## not updated yet ##########################
function LaplaceMMS(simul::SEM_Wave)   

    nx = length(simul.x)
    ny = length(simul.y)
    refLap = zeros(ny, nx) 
    ref = zeros(ny, nx) 

    # update the c_square matrix:
    for j = 1:length(simul.x)
        for i = 1:length(simul.y)
            simul.c_square[i, j] = 1.0 + MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0, 0, simul.MMS_c)
        end
    end

    



    for i = 1:ny
        for j = 1:nx
            ref[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 0, 0, simul.MMS_j)

            refLap[i, j] = simul.c_square[i, j] * MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 2, 0, 0, simul.MMS_j) + # Laplacian!
                           simul.c_square[i, j] * MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 2, 0, simul.MMS_j) + MMS.MMSfun(simul.x[j], simul.y[end+1-i], 1, 0, simul.MMS_c) * MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 1, 0, 0, simul.MMS_j) + 
                           MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0, 1, simul.MMS_c) * MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 1, 0, simul.MMS_j)

            #println(simul.c_square[i, j])
            #println("index: (" * string(i) * ", " * string(j) * "), point: (" * string(round(simul.x[j], sigdigits = 3)) * ", " * string(round(simul.y[end+1-i], sigdigits = 3)) * "), value: " * string(simul.c_square[i, j] * MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 2, 0, 0, simul.MMS_j)) * ", " * string(simul.c_square[i, j] * MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 2, 0, simul.MMS_j)) * " c: " * string(simul.c_square[i, j]) * " u_xx: " * string(MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 2, 0, 0, simul.MMS_j)) * " u_yy: " * string(MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 2, 0, simul.MMS_j)))
            #println("index: (" * string(i) * ", " * string(j) * "), point: (" * string(round(simul.x[j], sigdigits = 3)) * ", " * string(round(simul.y[end+1-i], sigdigits = 3)) * "), value: " * string(refLap[i, j]) * " c: " * string(simul.c_square[i, j]) * " u_xx: " * string(MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 2, 0, 0, simul.MMS_j)) * " u_yy: " * string(MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 2, 0, simul.MMS_j)))
            
            #println("index: (" * string(i) * ", " *string(j) * ") " * string(simul.c_square[i, j]))
            #println("index: (" * string(i) * ", " *string(j) * ") " * string(MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 2, 0, 0, simul.MMS_j)) * ", " * string(MMS.MMSfun(simul.x[j], simul.y[end+1-i], 0.0, 0, 2, 0, simul.MMS_j)))
        end
    end



    # construct the exact value of g for the MMS case
    gVals = zeros(ny, nx)    
    refNormalDer, refTimeDer = SEM_Wave_2d.BoundaryDersMMS(simul, 0.0)

    alpha = simul.bc[1]; beta = simul.bc[2]
    gVals = alpha * refTimeDer + beta .* simul.c_square .* refNormalDer

    #boundaryVals = simul.M_b .* refNormalDer
    boundaryVals = simul.M_b .* gVals

    laplaceVals = SEM_Wave_2d.LaplaceTerm(simul, ref) / simul.timestep^2 + boundaryVals ./ (simul.M)

    return laplaceVals, refLap

end

########################## not updated yet ##########################
function ForcingMMS(simul::SEM_Wave, stepnumber::Int64)   

    # only works if (simul.useMMS == true)

    t = (stepnumber-1)*simul.timestep
    fRef = zeros(length(simul.y), length(simul.x)) 

    for i = 1:length(simul.y)
        for j = 1:length(simul.x)
            fRef[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 2, 0, 0, simul.MMS_j) + MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 0, 2, 0, simul.MMS_j) - MMS.MMSfun(simul.x[j], simul.y[end+1-i], t, 0, 0, 2, simul.MMS_j)
        end
    end

    fVals = -SEM_Wave_2d.ForcingTerm(simul, stepnumber) / simul.timestep^2

    return [fVals, fRef]

end


########################## not updated yet ##########################
function BoundaryMMS(simul::SEM_Wave, stepnumber::Int64, bc::Vector{Float64})   

    # computes the line integrals of gMMS and g from BoundaryTerm along the boundary of our domain.

    # only works if (simul.useMMS == true)
    t = (stepnumber-1)*simul.timestep
    gRef = zeros(length(simul.y), length(simul.x)) 

    # if we use MMS, g should just be whatever the target solution has 
    refNormalDer = zeros(length(simul.y), length(simul.x))
    refTimeDer = zeros(length(simul.y), length(simul.x))
    
    simul.bc = bc

    alpha = simul.bc[1]; beta = simul.bc[2]

    for i = 1:length(simul.x)
        # upper side:
        refNormalDer[1, i] = MMS.MMSfun(simul.x[i], simul.y[end], t, 0, 1, 0, simul.MMS_j)
        refTimeDer[1, i] = MMS.MMSfun(simul.x[i], simul.y[end], t, 0, 0, 1, simul.MMS_j)

        # lower side: 
        refNormalDer[end, i] = - MMS.MMSfun(simul.x[i], simul.y[1], t, 0, 1, 0, simul.MMS_j)
        refTimeDer[end, i] = MMS.MMSfun(simul.x[i], simul.y[1], t, 0, 0, 1, simul.MMS_j)
    end

    for j = 1:length(simul.y)
        # left side: 
        refNormalDer[j, 1] = - MMS.MMSfun(simul.x[1], simul.y[end+1-j], t, 1, 0, 0, simul.MMS_j)
        refTimeDer[j, 1] = MMS.MMSfun(simul.x[1], simul.y[end+1-j], t, 0, 0, 1, simul.MMS_j)

        # right side: 
        refNormalDer[j, end] = MMS.MMSfun(simul.x[end], simul.y[end+1-j], t, 1, 0, 0, simul.MMS_j)
        refTimeDer[j, end] = MMS.MMSfun(simul.x[end], simul.y[end+1-j], t, 0, 0, 1, simul.MMS_j)
    end

    gRef = (alpha * refTimeDer + beta .* simul.c_square .* refNormalDer)

    # take data from BoundaryTerm
    gVals = (beta / simul.timestep^2) * (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b) .* SEM_Wave_2d.BoundaryTerm(simul, stepnumber) # not yet multiplied by simul.M_b

    gVals[:, 1] = gVals[:, 1] ./ simul.M_b[:, 1]
    gVals[:, end] = gVals[:, end] ./ simul.M_b[:, end]
    gVals[1, 2:end-1] = gVals[1, 2:end-1] ./ simul.M_b[1, 2:end-1]
    gVals[end, 2:end-1] = gVals[end, 2:end-1] ./ simul.M_b[end, 2:end-1]
    
    return [gVals, gRef]
    

end



function InitialiseMMS(simul::SEM_Wave, uStart::Matrix{Float64}, uStartDer::Matrix{Float64}, Tend::Float64, nsteps::Int64, forcing::Matrix{Float64}, omega::Float64, alphas::Vector{Float64}, g::Matrix{Float64})

    Initialise!(simul, uStart, uStartDer, Tend, nsteps, forcing, omega, alphas, g)

    uPrevRef = zeros(length(simul.y), length(simul.x))

    for j = 1:length(simul.x)
        for i = 1:length(simul.y)
            uPrevRef[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], -simul.timestep, 0, 0, 0, simul.MMS_j)
        end
    end

    return [simul.uPrev, uPrevRef]

end


function Waveholtz(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, alphas, g, tol::Float64)

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))
    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))



    ############################################
    #timestep = Tend/nsteps
    # now make that one more timestep, so that the central difference approximation of \partial_t u works at t = Tend
    #Tend = Tend + timestep
    #nsteps = nsteps + 1
    ############################################

    #println("Number of timesteps for waveholtz: " * string(nsteps) * " with a step size of " * string(timestep))

    # starting guess
    uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = zeros(length(simul.y), length(simul.x))

    animate = false     # we do not want any animations of the wave eq solutions
    res = Inf


    nIter = 0
    maxIter = 1e5


    while res > tol

        #to = TimerOutput()

        SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

        #show(to)

        res = SEM_Wave_2d_Updated.SeminormWH(simul, uStart - simul.uFiltered, uStartDer - simul.uDerFiltered) / SEM_Wave_2d_Updated.SeminormWH(simul, uStart, uStartDer)

        #println(maximum(abs.(simul.uFiltered)))
        #println(maximum(abs.(simul.uDerFiltered)))


        #res = (SEM_Wave_2d.GradIntegral(simul, uStart - simul.uFiltered, uStart - simul.uFiltered) + SEM_Wave_2d.LpNorm(simul, uStartDer, simul.uDerFiltered, 2)^2)^(1/2) / 
        #        (SEM_Wave_2d.GradIntegral(simul, simul.fVals, simul.fVals))^(1/2)



        #=
        # relative residual
        if nIter > 1

            res = (LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)^2 + 
                LpNorm(simul, simul.uDerFiltered, oldAppx[length(simul.y) + 1:end, :], 2)^2)^(0.5) / LpNorm(simul, simul.fVals, zeros(length(simul.y), length(simul.y)), 2)

        end
        =#

        uStart .= simul.uFiltered
        uStartDer .= simul.uDerFiltered
        

        nIter = nIter + 1

        if nIter > maxIter
            res = 0.0 # i.e. break the run
            println("maxIter reached for omega = " * string(omega))
        end

        #println("$omega || iteration: " * string(nIter) * " || residual: " * string(res))

    end


    return simul.uFiltered, simul.uDerFiltered, nIter

end


function WaveholtzGMRES(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, alphas, g, tol)

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))

    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))

    timestep = Tend/nsteps


    ############################################
    # now make that one more timestep, so that the central difference approximation of \partial_t u works at t = Tend
    #Tend = Tend + timestep
    #nsteps = nsteps + 1
    ############################################

    
    
    #println("$nsteps steps used for the GMRES-accelerated WaveHoltz")
    #println("time step size for the wave solver: " * string(timestep))

    nx = length(simul.x)
    ny = length(simul.y)
    N = Int(length(simul.x) * length(simul.y)) # half of the number of degrees of freedoms of our system
    

    # perfrom one WH step on the zero vector, this will be the RHS of our linear system
    uStart = zeros(ny, nx)
    uStartDer = zeros(ny, nx)

    SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)

    b = [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)]


    function WaveholtzAction(vec)    

        SEM_Wave_2d_Updated.Simulate(simul, Matrix(reshape(vec[1:Int(N)], ny, nx)), Matrix(reshape(vec[(Int(N)+1):end], ny, nx)), Tend, nsteps, fVals, omega, alphas, g, false)        
        w = vec - [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)] + b    

        return w
        
    end

    WaveholtzMatVec = LinearMap(WaveholtzAction, 2*N) # the matvec as a LinearMap

    #x, history = gmres(WaveholtzMatVec, b, verbose=true)
    #x, history = gmres(WaveholtzMatVec, b, log=true, verbose=false, reltol=tol,  restart=1000)    
    x, history = gmres(WaveholtzMatVec, b, log=true, verbose=true, reltol=tol,  restart=10000)    

    u_0 = Matrix(reshape(x[1:N, 1], ny, nx))
    u_1 = Matrix(reshape(x[N+1:end, 1], ny, nx))

    return u_0, u_1, history

end


function WaveholtzGMRES(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, alphas, g, nIter::Int64)

    # almost entirely a copy of the other WaveHoltzGMRES-function. Updates should be made in the near future

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))
    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))
    timestep = Tend/nsteps

    nx = length(simul.x)
    ny = length(simul.y)
    N = Int(length(simul.x) * length(simul.y)) # half of the number of degrees of freedoms of our system
    

    # perfrom one WH step on the zero vector, this will be the RHS of our linear system
    uStart = zeros(ny, nx)
    uStartDer = zeros(ny, nx)

    SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, false)

    b = [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)]


    function WaveholtzAction(vec)    

        SEM_Wave_2d_Updated.Simulate(simul, Matrix(reshape(vec[1:Int(N)], ny, nx)), Matrix(reshape(vec[(Int(N)+1):end], ny, nx)), Tend, nsteps, fVals, omega, alphas, g, false)        
        w = vec - [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)] + b    

        return w
        
    end

    WaveholtzMatVec = LinearMap(WaveholtzAction, 2*N) # the matvec as a LinearMap

    #x, history = gmres(WaveholtzMatVec, b, verbose=true)
    #x, history = gmres(WaveholtzMatVec, b, log=true, verbose=false, reltol=tol,  restart=1000)    
    x, history = gmres(WaveholtzMatVec, b, restart = nIter, maxiter = nIter, reltol = 0.0, abstol = 0.0, log = true, verbose = true)

    u_0 = Matrix(reshape(x[1:N, 1], ny, nx))
    u_1 = Matrix(reshape(x[N+1:end, 1], ny, nx))

    return u_0, u_1, history

end

function WaveholtzAnimation(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, alphas, g, maxIter::Int64, logscale = false)

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))
    nsteps = Integer(ceil(1.5*Tend * cMax * (1/delta_x + 1/delta_y)))
    timestep = Tend/nsteps

    println("time step size for the wave solver: " * string(timestep))

    # starting guess
    uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = zeros(length(simul.y), length(simul.x))
    oldAppx = zeros(2*length(simul.y), length(simul.x)) # stores uFiltered and uFilteredDer from the previous iteration

    animate = false     # we do not want any animations of the wave eq solutions
    res = NaN
    
    nIter = 1

    anim = Animation()

    
    
    if logscale
        surface(simul.x, simul.y[end:-1:1], log10.(abs.(simul.uFiltered)), legend=:false)
    else
        surface(simul.x, simul.y[end:-1:1], simul.uFiltered, zlims=(-1, 1), legend=:false)
    end
    frame(anim)


    while maxIter+1 > nIter
        
        println("iteration: " * string(nIter) * " || residual: " * string(res))
        
        SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered
        #uStartDer = zeros(length(simul.y), length(simul.x))


        # relative residual
        res = SEM_Wave_2d_Updated.SeminormWH(simul, uStart - simul.uFiltered, uStartDer - simul.uDerFiltered) / SEM_Wave_2d_Updated.SeminormWH(simul, uStart, uStartDer)

        #res = LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)/LpNorm(simul, simul.uFiltered, zeros(length(simul.y), length(simul.y)), 2)

        oldAppx = [simul.uFiltered; simul.uDerFiltered]

        nIter = nIter + 1

        if logscale
            plotHigh = maximum(log10.(abs.(simul.uFiltered)))
            plotLow = -15
            surface(simul.x, simul.y[end:-1:1], log10.(abs.(simul.uFiltered)), zlims=(plotLow, plotHigh), legend=:false)
        else
            plotHigh = maximum(simul.uFiltered)
            plotLow = minimum(simul.uFiltered)
            if plotHigh < 0.0
                plotHigh = 0.9*plotHigh
                plotLow = 1.1*plotLow
            elseif plotLow > 0.0
                plotHigh = 1.1*plotHigh
                plotLow = 0.9*plotLow
            else
                plotHigh = 1.1*plotHigh
                plotLow = 1.1*plotLow
            end
            
            #surface(simul.x, simul.y[end:-1:1], simul.uFiltered, zlims=(plotLow, plotHigh), legend=:false)
            surface(simul.x, simul.y[end:-1:1], simul.uFiltered, zlims=(-0.5, 0.5), legend=:false)
        end
        
        frame(anim) 

    end

    gif(anim, "WaveholtzGif.gif", fps=10)

end

function WaveholtzConvHistory(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, alphas, g, maxIter::Int64, makePlot::Bool = true)

    ### produces a plot of the convergence history of nIter Waveholtz iterations

    # plots the convergence history of a Waveholtz iteration
    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))

    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))

    timestep = Tend/nsteps

    println("Number of timesteps for waveholtz: " * string(nsteps) * " with a step size of " * string(timestep))
    

    # starting guess
    uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = zeros(length(simul.y), length(simul.x))
    oldAppx = zeros(2*length(simul.y), length(simul.x)) # stores uFiltered and uFilteredDer from the previous iteration

    animate = false     # we do not want any animations of the wave eq solutions
    res = NaN

    nIter = 1

    data = zeros(maxIter, 1)

    while maxIter+1 > nIter

        SEM_Wave_2d_Updated.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, alphas, g, animate)

        #fValsH1Norm = (SEM_Wave_2d.LpNorm(simul, simul.fVals, zeros(length(simul.y), length(simul.x)), 2)^2 + SEM_Wave_2d.GradIntegral(simul, simul.fVals, simul.fVals))^(1/2)
        #res = (SEM_Wave_2d.GradIntegral(simul, uStart - simul.uFiltered, uStart - simul.uFiltered) + SEM_Wave_2d.LpNorm(simul, uStartDer, simul.uDerFiltered, 2)^2)^(1/2) / fValsH1Norm


        res = SEM_Wave_2d_Updated.SeminormWH(simul, uStart - simul.uFiltered, uStartDer - simul.uDerFiltered) / SEM_Wave_2d_Updated.SeminormWH(simul, uStart, uStartDer)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

        # seminorm in the residual relative to the H1-norm of fVals
                

        #=
        res = (LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)^2 + 
               LpNorm(simul, simul.uDerFiltered, oldAppx[length(simul.y) + 1:end, :], 2)^2)^(0.5) /
              (LpNorm(simul, oldAppx[1:length(simul.y), :], zeros(length(simul.y), length(simul.x)), 2)^2 + 
               LpNorm(simul, oldAppx[length(simul.y) + 1:end, :], zeros(length(simul.y), length(simul.x)), 2)^2)^0.5
        =#

        data[nIter] = res

        #res = LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)/LpNorm(simul, simul.uFiltered, zeros(length(simul.y), length(simul.y)), 2)

        oldAppx = [simul.uFiltered; simul.uDerFiltered]

        nIter = nIter + 1

        println("iteration: " * string(nIter) * " || residual: " * string(res))

    end

    #plt = plot(collect(LinRange(1, maxIter, maxIter)), xscale=:log10, yscale=:log10, data)
    if makePlot
        plt = plot(collect(LinRange(1, maxIter, maxIter)), data, yscale=:log10)
        savefig(plt, "convHistory.png")
    end

    return oldAppx, data

end


function HelmholtzMatrix(simul::SEM_Wave, omega::Float64, alphas, fVals)
    # returns the matrix corresponding to the discrete Helmholtz operator (vectorised!)


    # setup boundary conditions and frequency
    SetBCmat!(simul, alphas)
    simul.omega = omega

    nx = length(simul.x)
    ny = length(simul.y)

    Kx = length(simul.xNodes) - 1
    Ky = length(simul.yNodes) - 1

    DoFs = Int(nx*ny)

    H = spzeros(DoFs, DoFs)

    a = simul.stepCoeffs[1]
    b = simul.stepCoeffs[2]

    e_j = zeros(DoFs, 1)
    Le_j = zeros(ny, nx)        


    for j = 1:DoFs

        prevIndex = maximum([j-1, 1])
        e_j[prevIndex] = 0.0
        e_j[j] = 1.0

        # appropriate shape for use with the calculations of the Laplacian
        e_j = reshape(e_j, ny, nx)
        
        # matrix for storing the data 
        Le_j .= zeros(size(e_j))        
        
        # coordinates of the nonzero element in the matrix e_j
        idx = findfirst(!iszero, e_j)
        row, col = idx.I

        # indices of elements we will compute the Laplacian in
        rowIndexList = [] 
        colIndexList = [] 

        if (mod(row, simul.N) == 1) # if the nonzero point is on the border of two elements in the y-direction

            if (row == 1)

                rowIndex = Int((row-1 - mod(row-1, simul.N))/simul.N)+1
                rowIndexList = [rowIndexList; rowIndex]

                
            elseif (row == ny) # if on edge of physical domain

                rowIndex = Int((row-1 - mod(row-1, simul.N))/simul.N)
                rowIndexList = [rowIndexList; rowIndex]
                
            else # there must be two relevant row indices, namely

                rowIndex1 = Int((row-1 - mod(row-1, simul.N))/simul.N) + 1
                rowIndex2 = rowIndex1 - 1

                rowIndexList = [rowIndexList; rowIndex1; rowIndex2]

            end

        else

            rowIndex = Int((row - 1 - mod(row-1, simul.N))/simul.N) + 1
            rowIndexList = [rowIndexList; rowIndex]

        end


        if (mod(col, simul.N) == 1) # if the nonzero point is on the border of two elements in the x-direction

            if (col == 1)

                colIndex = Int((col-1 - mod(col-1, simul.N))/simul.N) + 1
                colIndexList = [colIndexList; colIndex]


            elseif (col == nx) # if on edge of physical domain

                colIndex = Int((col-1 - mod(col-1, simul.N))/simul.N)
                colIndexList = [colIndexList; colIndex]


            else # there must be two relevant column indices, namely

                colIndex1 = Int((col-1 - mod(col-1, simul.N))/simul.N) + 1
                colIndex2 = colIndex1 - 1

                colIndexList = [colIndexList; colIndex1; colIndex2]

            end

        else
            colIndex = Int((col-1 - mod(col-1, simul.N))/simul.N) + 1
            colIndexList = [colIndexList; colIndex]
        end


        for s = 1:length(rowIndexList) # at most 4 relevant elements for e_j, we consider each one.
            for t = 1:length(colIndexList)

                rowIndex = rowIndexList[s]
                colIndex = colIndexList[t]

                # index appropriate for the LocalLaplace function
                #k = ((Ky - rowIndex) * Kx) + colIndex # my idea
                k = ((rowIndex - 1) * Kx) + colIndex   # Claude fix

                #k = Int(k)


                # perform the Laplace calculation
                SEM_Wave_2d_Updated.LocalLaplace!(simul, k, e_j, Le_j)

            end
        end


        # store data
        H[:, j] .= reshape(Le_j, DoFs, 1)

        println("generated row $j out of $DoFs")

    end

    M = spdiagm(reshape(simul.M, DoFs, 1)[:])

    B = spdiagm(reshape(simul.bcMat .* simul.M_b, DoFs, 1)[:])

    H = simul.omega^2 * M + (1im*simul.omega)*B + H

    F = M * reshape(fVals, DoFs, 1)

    return H, F

end


function HelmholtzMatrixOld(simul::SEM_Wave, omega::Float64, alphas, fVals)
    # returns the matrix corresponding to the discrete Helmholtz operator (vectorised!)


    # setup boundary conditions and frequency
    SetBCmat!(simul, alphas)
    simul.omega = omega

    nx = length(simul.x)
    ny = length(simul.y)

    DoFs = Int(nx*ny)

    H = spzeros(DoFs, DoFs)
    
    a = simul.stepCoeffs[1]
    b = simul.stepCoeffs[2]


    for j = 1:DoFs

        e_j = zeros(DoFs, 1)
        e_j[j] = 1.0

        # This solution is extremely inefficient, and gets relatively worse the larger the problem is. 
        # Almost every element is full of zeros, so most of the work done in LaplaceTerm is adding and multiplying zeros. 
        Le_j = sparse(((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* SEM_Wave_2d_Updated.LaplaceTerm(simul, reshape(e_j, ny, nx)) / simul.timestep^2)

        # store data
        H[:, j] = reshape(Le_j, DoFs, 1)

        println("generated row $j out of $DoFs")

    end

    M = spdiagm(reshape(simul.M, DoFs, 1)[:])
    
    B = spdiagm(reshape(simul.bcMat .* simul.M_b, DoFs, 1)[:])

    H = simul.omega^2 * M + (1im*simul.omega)*B + H

    F = M * reshape(fVals, DoFs, 1)

    return H, F

end


function LocalLaplace!(simul::SEM_Wave, k::Int64, U::Matrix{Float64}, laplaceVals::Matrix{Float64})
    # approximates the Laplacian of U in element k and stores the result in laplaceVals
    # intended to replace the corresponding part of the code in LaplaceTerm by a function call

    pointsPerElement = simul.N + 1
    V_k = zeros(pointsPerElement, pointsPerElement)
    U_k = GetDegreesOfFreedom(simul, k, U)
    #c_square_k = GetDegreesOfFreedom(simul, k, simul.c_square) # no longer needed

    V_kx = zeros(pointsPerElement, pointsPerElement)
    V_ky = zeros(pointsPerElement, pointsPerElement)

    for j = 1:pointsPerElement # fill the matrices, column by column or row by row:

        mul!((@view V_kx[j, :]), simul.DxMatrices[(k-1)*pointsPerElement + j], @view U_k[j, :])
        mul!((@view V_ky[:, j]), simul.DyMatrices[(k-1)*pointsPerElement + j], @view U_k[:, j])

    end


    V_kx .= simul.QuadWeights .* V_kx
    V_ky .= simul.QuadWeights' .* V_ky
    V_k .= V_kx .+ V_ky

    V_k .= -V_k

    # store the result in laplaceVals
    SetDegreesOfFreedom!(simul, k, laplaceVals, V_k, true)

end


# does not work as intended yet!
function DirectSolveGMRES(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, alphas, g, tol)


    println("the function DirectSolveGMRES does not work as intended!")
    nx = length(simul.x)
    ny = length(simul.y)
    
    simul.omega = omega
    simul.fVals = fVals
    simul.alphas = alphas
    
    SetBCmat!(simul, alphas)

    simul.g = g

    a = simul.stepCoeffs[1]
    b = simul.stepCoeffs[2]

    nx = length(simul.x)
    ny = length(simul.y)
    N = Int(nx*ny)

    function HelmholtzAction(vec)    

        U = Matrix(reshape(vec, ny, nx))
        
        HU = omega^2 .* (simul.M .* U) + 
                        (1im * omega) .* ((simul.bcMat) .* simul.M_b .* U) + 
                        ((simul.M ./ a) + (simul.timestep ./ (2*b)) * simul.bcMat .* simul.M_b) .* SEM_Wave_2d_Updated.LaplaceTerm(simul, U) / simul.timestep^2
        

        #println(typeof((1im * omega) .* ((simul.bcMat) .* simul.M_b .* U)))
        #println(maximum(abs.(imag(((1im * omega) .* ((simul.bcMat) .* simul.M_b .* U))))))
        
        #println(maximum(abs.(imag(HU))))

        return reshape(HU, N)
        
    end

    HelmholtzMatVec = LinearMap{ComplexF64}(HelmholtzAction, N) # the matvec as a linear map
    f = reshape(fVals, N)                           # right-hand side

    


    x0 = zeros(ComplexF64, N) # complex initial guess so that gmres returns complex solutions

    x, history = gmres(HelmholtzMatVec, f, log=true, verbose=true, reltol=tol,  restart=1000)    
    #x, history = gmres(HelmholtzMatVec, f, log=true, verbose=true, reltol=tol,  restart=Int(round(N/2)))    

    U_sol = Matrix(reshape(x, ny, nx))

    return U_sol, history

end




function ErrorEstimate(simul::SEM_Wave, u::Matrix{Float64}, fVals::Matrix{Float64}, omega::Float64) 
    
    ### should estimate || Laplace u + omega^2 u - f ||_{interior nodes} + || alpha u + \beta \omega \partial_n u ||_{boundary}, which should show whether u is a good approximation to the corresponding PDE solution

    println("ErrorEstimate is under construction: do not use it unless you understand the details!")

    laplaceVals = SEM_Wave_2d.LaplaceTerm(simul, u) / simul.timestep^2 # works since we will not look at the boundary, meaning that we can ignore the boundary integral

    #error_1 = LpNorm(simul, laplaceVals[2:end-1, 2:end-1], fVals[2:end-1, 2:end-1]) this might actually add some error corresponding to the boundary bits
    error_1 = norm(laplaceVals[2:end-1, 2:end-1] .+ omega^2 .* u[2:end-1, 2:end-1] .- fVals[2:end-1, 2:end-1], 2) ./ norm(fVals[2:end-1, 2:end-1], 2)
    #println(norm(laplaceVals[2:end-1, 2:end-1] .+ omega^2 .* u[2:end-1, 2:end-1] .- fVals[2:end-1, 2:end-1], 2))

    return error_1

end



function single_gs(Q, w, k)
    
    h = Q[:, 1:k]'*w;
    
    w_orth = w - Q[:, 1:k]*h;
    
    beta = norm(w_orth);

    return h, beta, w_orth

end


function double_gs(Q, w, k)

    h, beta, w_orth = single_gs(Q, w, k);
    
    g, beta, w_orth = single_gs(Q, w_orth, k);
    
    h = h + g;

    return h, beta, w_orth

end


function LpNorm(simul::SEM_Wave, u::Matrix{Float64}, v::Matrix{Float64}, p::Int64)

    ### calculates the Lp norm of u - v using LGL quadrature; u and v have values for each point given by simul.x, simul.y

    pointsPerElement = simul.N + 1                                 # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        
    
   
    QuadWeightsMat = simul.QuadWeights * simul.QuadWeights' 

    integrand_k = zeros(pointsPerElement, pointsPerElement)

    norm = 0.0

    for k = 1:Kx*Ky

        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = simul.xNodes[i+1] - simul.xNodes[i]
        delta_y_k = simul.yNodes[j+1] - simul.yNodes[j]

        integrand_k .= (abs.(SEM_Wave_2d_Updated.GetDegreesOfFreedom(simul, k, u) .- SEM_Wave_2d_Updated.GetDegreesOfFreedom(simul, k, v))).^p

        
        # integrate
        norm = norm + 0.25*(delta_x_k*delta_y_k) * sum(QuadWeightsMat .* integrand_k)

    end

    norm = norm^(1/p)

    return norm

end

function SeminormWH(simul::SEM_Wave, U_real::Matrix{Float64}, U_imag::Matrix{Float64})

    # calculates that seminorm that the convergence of Waveholtz is monotone in

    #println(SEM_Wave_2d.GradIntegral(simul, U_real, U_real)) # sometimes this becomes a small negative number, especially when the polynomial order is high.
    #println(SEM_Wave_2d.LpNorm(simul, U_imag, zeros(length(simul.y), length(simul.x)), 2))

    #the abs is there for stability reasons: sometimes the GradIntegral becomes -1e-14 or similar, and breaks this square root. I think this is due to float errors, as it gets worse with higher orders...
    out = (abs(SEM_Wave_2d_Updated.GradIntegralSeminorm(simul, U_real, U_real) + SEM_Wave_2d_Updated.LpNorm(simul, U_imag, zeros(length(simul.y), length(simul.x)), 2)^2))^(0.5)

    return out

end


function GradIntegralSeminorm(simul::SEM_Wave, U::Matrix{Float64}, V::Matrix{Float64})
    ### calculates \int_\Omega c^2(x, y) \grad u \cdot \grad v dxdy where the nodal values of u and v are stored in the matrices U and V
    ### for use in SeminormWH


    pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        

    gradIntegral = 0.0                                                  # the final output will be a real number, with a contribution from each element

    U_k = zeros(pointsPerElement, pointsPerElement)                     # the values of U in element k
    V_k = zeros(pointsPerElement, pointsPerElement)                     # the values of V in element k
    G = simul.G
    D = simul.D

    W = diagm(simul.QuadWeights)

    for k = 1:Kx*Ky                                                     # loop over the elements and calculate \int_{\Omega_k} c^2(x, y) \grad u \cdot \grad v dxdy

        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = simul.xNodes[i+1] - simul.xNodes[i]
        delta_y_k = simul.yNodes[j+1] - simul.yNodes[j]

        U_k .= GetDegreesOfFreedom(simul, k, U)
        V_k .= GetDegreesOfFreedom(simul, k, V)
        c_square_k = GetDegreesOfFreedom(simul, k, simul.c_square) # superfluous?

        xDerContribution = 0.0
        yDerContribution = 0.0

        for t = 1:pointsPerElement
            
            wt = simul.QuadWeights[t]

            #At = transpose(simul.D) * W * diagm(c_square_k[t, :]) * simul.D
            At = simul.DxMatrices[(k-1)*pointsPerElement + t]
            #Bt = transpose(simul.D) * W * diagm(c_square_k[:, t]) * simul.D
            Bt = simul.DyMatrices[(k-1)*pointsPerElement + t]

            xSum = 0.0
            ySum = 0.0

            for l = 1:pointsPerElement
                for m = 1:pointsPerElement 
                    xSum = xSum + wt * U_k[l, t] * V_k[m, t] * At[l, m]
                    ySum = ySum + wt * U_k[t, l] * V_k[t, m] * Bt[l, m]
                end
            end

            xDerContribution = xDerContribution + (delta_y_k/delta_x_k) * xSum
            yDerContribution = yDerContribution + (delta_x_k/delta_y_k) * ySum

        end

        gradIntegral_k = xDerContribution + yDerContribution
        
        gradIntegral = gradIntegral + gradIntegral_k

    end


    # I have had situations where this returns negative numbers of the order 1e-17
    # I suppose these are zero up to flop-errors, and set them to zero in that case, so that Seminorm and SobolevNorm do not have issues down the line
    if abs(gradIntegral) < 1e-12 # importantly, if there were a bug that made gradIntegral < -1e-12, say, this would not hide that problem.
        gradIntegral = 0.0
    end

    return gradIntegral

end

function GradIntegral(simul::SEM_Wave, U::Matrix{Float64}, V::Matrix{Float64})
    ### calculates \int_\Omega \grad u \cdot \grad v dxdy where the nodal values of u and v are stored in the matrices U and V
    ### for use in SobolevNorm


    pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        

    gradIntegral = 0.0                                                  # the final output will be a real number, with a contribution from each element

    U_k = zeros(pointsPerElement, pointsPerElement)                     # the values of U in element k
    V_k = zeros(pointsPerElement, pointsPerElement)                     # the values of V in element k
    G = simul.G
    D = simul.D
    

    W = diagm(simul.QuadWeights)

    for k = 1:Kx*Ky                                                     # loop over the elements and calculate \int_{\Omega_k} \grad u \cdot \grad v dxdy


        c_square_k = GetDegreesOfFreedom(simul, k, simul.c_square) # superfluous?

        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = simul.xNodes[i+1] - simul.xNodes[i]
        delta_y_k = simul.yNodes[j+1] - simul.yNodes[j]

        U_k .= GetDegreesOfFreedom(simul, k, U)
        V_k .= GetDegreesOfFreedom(simul, k, V)

        xDerContribution = 0.0
        yDerContribution = 0.0

        for t = 1:pointsPerElement

            wt = simul.QuadWeights[t]


            # like the matrices in simul.DxMatrices, simul.DyMatrices, but with c^2 = 1.0...
            At = (delta_y_k/delta_x_k) * transpose(simul.D) * W * simul.D
            Bt = (delta_x_k/delta_y_k) * transpose(simul.D) * W * simul.D
            
            xSum = 0.0
            ySum = 0.0

            for l = 1:pointsPerElement
                for m = 1:pointsPerElement 
                    xSum = xSum + wt * U_k[l, t] * V_k[m, t] * At[l, m]
                    ySum = ySum + wt * U_k[t, l] * V_k[t, m] * Bt[l, m]
                end
            end

            xDerContribution = xDerContribution + (delta_y_k/delta_x_k) * xSum
            yDerContribution = yDerContribution + (delta_x_k/delta_y_k) * ySum

        end

        gradIntegral_k = xDerContribution + yDerContribution
        gradIntegral = gradIntegral + gradIntegral_k

    end


    # I have had situations where this returns negative numbers of the order 1e-17
    # I suppose these are zero up to flop-errors, and set them to zero in that case, so that Seminorm and SobolevNorm do not have issues down the line
    if abs(gradIntegral) < 1e-16 # importantly, if there were a bug that made gradIntegral < -1e-12, say, this would not hide that problem.
        gradIntegral = 0.0
    end

    return gradIntegral

end

function SobolevNormOld(simul::SEM_Wave, U::Matrix{Float64}, k::Float64)

    ### calculates the H_1^k-norm of u whose values is stored in the matrix U, i.e. \|U\|_{L^2} + k^{-2}\|\nabla U\|_{L^2}

    pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        

    gradIntegral = 0.0                                                  # the final output will be a real number, with a contribution from each element

    U_k = zeros(pointsPerElement, pointsPerElement)                     # the values of U in element k
    G = simul.G
    D = simul.D

    W = diagm(simul.QuadWeights)

    for k = 1:Kx*Ky                                                     # loop over the elements and calculate \int_{\Omega_k}  |\grad u|^2 dxdy

        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = simul.xNodes[i+1] - simul.xNodes[i]
        delta_y_k = simul.yNodes[j+1] - simul.yNodes[j]

        U_k .= GetDegreesOfFreedom(simul, k, U)

        xDerContribution = 0.0
        yDerContribution = 0.0

        for t = 1:pointsPerElement
            
            wt = simul.QuadWeights[t]
            DtWD = transpose(simul.D) * W * simul.D

            xSum = 0.0
            ySum = 0.0

            for l = 1:pointsPerElement
                for m = 1:pointsPerElement 
                    xSum = xSum + wt * U_k[l, t] * U_k[m, t] * DtWD[l, m]
                    ySum = ySum + wt * U_k[t, l] * U_k[t, m] * DtWD[l, m]
                end
            end

            xDerContribution = xDerContribution + (delta_y_k/delta_x_k) * xSum
            yDerContribution = yDerContribution + (delta_x_k/delta_y_k) * ySum

        end

        gradIntegral_k = xDerContribution + yDerContribution
        
        gradIntegral = gradIntegral + gradIntegral_k

    end

    return LpNorm(simul, U, zeros(length(simul.y), length(simul.x)), 2) + sqrt(gradIntegral) / k^2


end

function SobolevNorm(simul::SEM_Wave, U::Matrix{Float64}, k::Float64)

    ### calculates the H_1^k-norm of u whose values is stored in the matrix U, i.e. \|U\|_{L^2} + k^{-2}\|\nabla U\|_{L^2}

    return (LpNorm(simul, U, zeros(length(simul.y), length(simul.x)), 2)^2 + GradIntegral(simul, U, U) / k^2)^(0.5)


end



function BoundaryDersMMS(simul::SEM_Wave, t::Float64)
    ### calculates the weighted normal and time derivatives of the MMS solution on the boundary at time t, for use in MMS tests
    ### by "weighted normal derivative" we here mean the quantity \beta * c * \partial_n

    #weightedRefNormalDer = zeros(length(simul.y), length(simul.x))
    refNormalDer = zeros(length(simul.y), length(simul.x))
    refTimeDer = zeros(length(simul.y), length(simul.x))
    betas = sqrt.(1 .- simul.alphas.^2)

    for i = 1:length(simul.x)
        # upper side:
        #weightedRefNormalDer[1, i] = betas[2] * simul.c_square[1, i] * MMS.MMSfun(simul.x[i], simul.y[end], t, 0, 1, 0, simul.MMS_j)
        refNormalDer[1, i] = MMS.MMSfun(simul.x[i], simul.y[end], t, 0, 1, 0, simul.MMS_j)
        refTimeDer[1, i] = MMS.MMSfun(simul.x[i], simul.y[end], t, 0, 0, 1, simul.MMS_j)

        # lower side: 
        #weightedRefNormalDer[end, i] = - betas[4] * simul.c_square[end, i] * MMS.MMSfun(simul.x[i], simul.y[1], t, 0, 1, 0, simul.MMS_j)
        refNormalDer[end, i] = - MMS.MMSfun(simul.x[i], simul.y[1], t, 0, 1, 0, simul.MMS_j)
        refTimeDer[end, i] = MMS.MMSfun(simul.x[i], simul.y[1], t, 0, 0, 1, simul.MMS_j)
    end

    for j = 1:length(simul.y)
        # left side: 
        #weightedRefNormalDer[j, 1] = - betas[3] * simul.c_square[j, 1] * MMS.MMSfun(simul.x[1], simul.y[end+1-j], t, 1, 0, 0, simul.MMS_j)
        refNormalDer[j, 1] = - MMS.MMSfun(simul.x[1], simul.y[end+1-j], t, 1, 0, 0, simul.MMS_j)
        refTimeDer[j, 1] = MMS.MMSfun(simul.x[1], simul.y[end+1-j], t, 0, 0, 1, simul.MMS_j)

        # right side: 
        #weightedRefNormalDer[j, end] = betas[1] * simul.c_square[j, end] * MMS.MMSfun(simul.x[end], simul.y[end+1-j], t, 1, 0, 0, simul.MMS_j)
        refNormalDer[j, end] = MMS.MMSfun(simul.x[end], simul.y[end+1-j], t, 1, 0, 0, simul.MMS_j)
        refTimeDer[j, end] = MMS.MMSfun(simul.x[end], simul.y[end+1-j], t, 0, 0, 1, simul.MMS_j)
    end


    # we need to take a weighted average in order to take account for the cusp at the corners
    delta_x = simul.xNodes[2] - simul.xNodes[1]
    delta_y = simul.yNodes[2] - simul.yNodes[1]

    #weightedRefNormalDer[1, 1] = simul.c_square[1, 1] * (-betas[3] * MMS.MMSfun(simul.x[1], simul.y[end], t, 1, 0, 0, simul.MMS_j)*delta_y + betas[2] * MMS.MMSfun(simul.x[1], simul.y[end], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) #upper left corner
    #weightedRefNormalDer[1, end] = simul.c_square[1, end] * (betas[1] * MMS.MMSfun(simul.x[end], simul.y[end], t, 1, 0, 0, simul.MMS_j)*delta_y + betas[2] * MMS.MMSfun(simul.x[end], simul.y[end], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y)  # upper right corner
    #weightedRefNormalDer[end, 1] = simul.c_square[end, 1] * (-betas[3] * MMS.MMSfun(simul.x[1], simul.y[1], t, 1, 0, 0, simul.MMS_j)*delta_y - betas[4] * MMS.MMSfun(simul.x[1], simul.y[1], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) # lower left corner
    #weightedRefNormalDer[end, end] = simul.c_square[end, end] * (betas[1] * MMS.MMSfun(simul.x[end], simul.y[1], t, 1, 0, 0, simul.MMS_j)*delta_y - betas[4] * MMS.MMSfun(simul.x[end], simul.y[1], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) # lower right corner

    refNormalDer[1, 1] = (-MMS.MMSfun(simul.x[1], simul.y[end], t, 1, 0, 0, simul.MMS_j)*delta_y + MMS.MMSfun(simul.x[1], simul.y[end], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 
    refNormalDer[1, end] = (MMS.MMSfun(simul.x[end], simul.y[end], t, 1, 0, 0, simul.MMS_j)*delta_y + MMS.MMSfun(simul.x[end], simul.y[end], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 
    refNormalDer[end, 1] = (-MMS.MMSfun(simul.x[1], simul.y[1], t, 1, 0, 0, simul.MMS_j)*delta_y - MMS.MMSfun(simul.x[1], simul.y[1], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 
    refNormalDer[end, end] = (MMS.MMSfun(simul.x[end], simul.y[1], t, 1, 0, 0, simul.MMS_j)*delta_y - MMS.MMSfun(simul.x[end], simul.y[1], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 

    
    #=
    println(weightedRefNormalDer[1, 1])
    println(weightedRefNormalDer[1, end])
    println(weightedRefNormalDer[end, 1])
    println(weightedRefNormalDer[end, end])
    =#

    return refNormalDer, refTimeDer
    #return weightedRefNormalDer, refTimeDer

end



end