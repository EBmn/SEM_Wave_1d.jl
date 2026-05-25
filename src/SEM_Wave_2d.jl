

module SEM_Wave_2d

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

    xNodes::Vector{Float64}         # the endpoints of the elements in the x-direction
    yNodes::Vector{Float64}         # the endpoints of the elements in the y-direction
    N::Int64                        # the degree of the Lagrange polynomials 
    x::Vector{Float64}              # the actual coordinates between xNodes[1] and xNodes[end] where we approximate our solution
    y::Vector{Float64}              # the actual coordinates between yNodes[1] and yNodes[end] where we approximate our solution
    c_square::Matrix{Float64}       # matrix storing the value of the squared wave speed at each point in simul.x * simul.y'

    nsteps::Int64                   # number of steps we make in time
    Tend::Float64                   # the endpoint of our interval in time
    timestep::Float64               # size of the timesteps (= Tend/nsteps)
    
    uPrev::Matrix{Float64}          # value of the displacement in the previous timepoint
    uNow::Matrix{Float64}           # value of the displacement in the current timepoint
    uNext::Matrix{Float64}          # value of the displacement in the next timepoint
    scratch::Matrix{Float64}        # scratch space
    uFiltered::Matrix{Float64}      # the filtered time average of u for use in Waveholtz
    uDerFiltered::Matrix{Float64}   # the filtered time average of u_t for use in Waveholtz

    fVals::Matrix{Float64}          # values of the driving term of the equation
    omega::Float64                  # frequency parameter omega
    bc::Vector{Float64}             # contains the \alpha and \beta used in the bc 
    g::Matrix{Float64}              # function values for g in Boundary condition
    
    
    QuadPoints::Vector{Float64}     # the points in [-1, 1] we use in our Gauss-Lobatto quadrature
    QuadWeights::Vector{Float64}    # the weights corresponding to the points in QuadPoints
    D::Array{Float64}               # matrix storing derivatives of Lagrange polynomials: D_{ij} = l'_j(\xi_i)
    G::Array{Float64}               # the reference stiffness matrix: G_{jm} = \sum_{l=0}^{N} w_l l'_m(\xi_l) l'_j(\xi_l); w_l \in QuadWeights, \xi_l \in QuadPoints, l_j's are the Lagrange polynomials                            
    M::Matrix{Float64}              # the mass matrix
    M_b::Matrix{Float64}            # the mass matrix on the boundary, used to handle terms from the boundary intergral
    #DxMatrices::Vector{SMatrix}
    #DyMatrices::Vector{SMatrix}
    DxMatrices::Vector{Matrix{Float64}}
    DyMatrices::Vector{Matrix{Float64}}
    
    useMMS::Bool                    # determines whether an MMS test is run or not
    MMS_j::MMS_jet                  # stores the relevant data for the MMS test
    MMS_c::MMS_jet                  # stores the relevant data for the MMS test with variable wave speed
    MMS_intDiff::Matrix{Float64}    # the integral of the absolute difference between the approximation and the MMS-reference at each point
    MMS_int::Matrix{Float64}        # the integral of the absolute value of the solution in time, for use for computing relative errors


    function SEM_Wave(
        xLims::Vector{Float64},     # smallest and largest x-values
        yLims::Vector{Float64},     # smallest and largest y-values
        Kx::Int64,                  # number of elements in x-direction
        Ky::Int64,                  # number of elements in y-direction
        N::Int64,                   # degree of polynomials; every element has N + 1 points
        c_square::Function          # wave speed
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
        bc = [0.0, 1.0]
        g = zeros(Ky*N + 1, Kx*N + 1)

        uFiltered = zeros(Ky*N + 1, Kx*N + 1)
        uDerFiltered = zeros(Ky*N + 1, Kx*N + 1)
        

        useMMS = false
        MMS_intDiff = zeros(Ky*N + 1, Kx*N + 1)
        MMS_int = zeros(Ky*N + 1, Kx*N + 1)

        Tend = 1                # arbitrary values, just here as placeholders
        nsteps = 100
        timestep = Tend/nsteps

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
                uPrev, uNow, uNext, scratch, uFiltered, uDerFiltered, fVals, omega, bc, g, QuadPoints, QuadWeights, D, G, M, M_b, DxMatrices, DyMatrices, useMMS, MMS_j, MMS_c, MMS_intDiff, MMS_int)
        
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
                  #LaplaceTermOld(simul, simul.uNow) + 
                  ForcingTerm(simul, stepnumber) + 
                  BoundaryTerm(simul, stepnumber)

    lap1 = LaplaceTermOld(simul, simul.uNow)
    lap2 = LaplaceTerm(simul, simul.uNow)
    println(maximum(abs.(lap1 - lap2)))

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

    ########################### OLD VERSION ###########################
    #=
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

    =#
    ###################################################################


    ########################### NEW VERSION ###########################
    
    @timeit to "time derivative" uDer = (simul.uNext - simul.uPrev) ./ (2*simul.timestep) # derivative at time t - simul.timestep (second order)

    # filter the solution for use in Waveholtz
    K = (cos(simul.omega*(t-simul.timestep)) - 0.25) * (2/simul.Tend)

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
        
        #println("should be (3*omega)/(4*pi) = " *string((3*simul.omega)/(4*pi)))
        #println(K)

        @timeit to "trapezoidal rule real part" simul.uFiltered = simul.uFiltered + simul.uNext * K * 0.5 * simul.timestep # simul.uNext corresponds to time t = simul.Tend
        @timeit to "trapezoidal rule imag part" simul.uDerFiltered = simul.uDerFiltered + uDer * K * 0.5 * simul.timestep 

    end
    
    

    # update the variables
    simul.uPrev = simul.uNow
    simul.uNow = simul.uNext

    #println(maximum(abs.(simul.uNow)))

    ###################################################################

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

function LaplaceTermOld(simul::SEM_Wave, U::Matrix{Float64})

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

    alpha = simul.bc[1]; beta = simul.bc[2]


    @timeit to "make" scratchMat = zeros(pointsPerElement, pointsPerElement)

    #@timeit to "make" V_kx = MMatrix{pointsPerElement, pointsPerElement, Float64}(undef)  # Mutable static matrix
    #@timeit to "make" V_ky = MMatrix{pointsPerElement, pointsPerElement, Float64}(undef)

    #=
    @timeit to "new make" row1 = zeros(pointsPerElement, 1)
    @timeit to "new make" row2 = zeros(pointsPerElement, 1)
    @timeit to "new make" col1 = zeros(pointsPerElement, 1)
    @timeit to "new make" col2 = zeros(pointsPerElement, 1)
    @timeit to "new make" scratch = zeros(pointsPerElement, 1)
    @timeit to "new make" scratchMat = zeros(pointsPerElement, pointsPerElement)
    =#

    for k = 1:Kx*Ky

        #=

        i = mod(k - 1, Kx) + 1
        j = Int((k - i) / Kx) + 1

        delta_x_k = simul.xNodes[i+1] - simul.xNodes[i]
        delta_y_k = simul.yNodes[j+1] - simul.yNodes[j]

        =#
       
        @timeit to "get" U_k .= GetDegreesOfFreedom(simul, k, U)

        @timeit to "get" c_square_k = GetDegreesOfFreedom(simul, k, simul.c_square)

        #@turbo for j = 1:pointsPerElement # fill the matrices, column by column or row by row:
        for j = 1:pointsPerElement # fill the matrices, column by column or row by row:

            ########### optimised(?) version ###########

            #=

            @timeit to "views" scratch = @view U_k[j, :]
            @timeit to "views" scratchMat = @views simul.D
            @timeit to "row*full" transpose(mul!(row1, scratchMat, scratch, delta_y_k/delta_x_k, 0)) # what is even the point of the transposition here?
            #@timeit to "row*full" transpose(mul!(row1, simul.D, scratch, delta_y_k/delta_x_k, 0))
            
            #@timeit to "row*full" transpose(mul!(row1, simul.D, U_k[j, :], delta_y_k/delta_x_k, 0))
            @timeit to "transposing" row1 = transpose(row1)
            #@timeit to "new row*diags" rmul!(row1, Diagonal(simul.QuadWeights .* c_square_k[j, :]))
            
            #=
            for t = 1:pointsPerElement
                @timeit to "row*diags" row1[t] = row1[t]*simul.QuadWeights[t]*c_square_k[j, t]
            end
            =#

            @timeit to "row*diags" rmul!(row1, Diagonal(simul.QuadWeights))
            @timeit to "row*diags" rmul!(row1, Diagonal(@view c_square_k[j, :]))

            #@timeit to "views" scratch = @views transpose(row1)
            @timeit to "views" row1 = transpose(@views row1)
            @timeit to "views" scratchMat = transpose(@views simul.D)
            #@timeit to "row*full" transpose(mul!(row2, scratchMat, row1))
            @timeit to "row*full" transpose(mul!(row2, scratchMat, scratch))
            #@timeit to "row*full" transpose(mul!(row2, transpose(simul.D), row1))
            @timeit to "store" V_kx[j, :] = row2
            
            @timeit to "views" scratch = @view U_k[:, j]
            @timeit to "full*col" transpose(mul!(col1, simul.D, scratch, delta_x_k/delta_y_k, 0))
            #@timeit to "new col*diags" lmul!(Diagonal(simul.QuadWeights .* c_square_k[:, j]), col1)

            @timeit to "col*diags" lmul!(Diagonal(simul.QuadWeights), col1)
            @timeit to "col*diags" lmul!(Diagonal(@view c_square_k[:, j]), col1)

            @timeit to "full*col" mul!(col2, scratchMat, col1)
            #@timeit to "full*col" mul!(col2, transpose(simul.D), col1)
            @timeit to "store" V_ky[:, j] = col2

            #@timeit to "xDer matmuls" V_kx[j, :] = U_k[j, :]' * (transpose(simul.D) * W * diagm(c_square_k[j, :]) * simul.D)
            #@timeit to "yDer matmuls" V_ky[:, j] = transpose(simul.D) * W * diagm(c_square_k[:, j]) * simul.D * U_k[:, j]

            =#


            ########### The Version to End All Versions ###########
            @timeit to "VecMat" mul!((@view V_kx[j, :]), simul.DxMatrices[(k-1)*pointsPerElement + j], @view U_k[j, :])
            @timeit to "MatVec" mul!((@view V_ky[:, j]), simul.DyMatrices[(k-1)*pointsPerElement + j], @view U_k[:, j])

            #index = (k-1)*pointsPerElement + j
            #@timeit to "VecMat2" scratchMat .= simul.DxMatrices[index]
            #@timeit to "VecMat" BLAS.gemv!('T', 1.0, scratchMat, @view(U_k[j, :]), 0.0, @view(V_kx[j, :]))
            #@timeit to "VecMat" BLAS.gemv!('T', 1.0, simul.DxMatrices[index], @view(U_k[j, :]), 0.0, @view(V_kx[j, :]))
            #@timeit to "MatVec2" scratchMat .= simul.DyMatrices[index]
            #@timeit to "MatVec" BLAS.gemv!('N', 1.0, scratchMat, @view(U_k[:, j]), 0.0, @view(V_ky[:, j]))
            #@timeit to "MatVec" BLAS.gemv!('N', 1.0, simul.DyMatrices[index], @view(U_k[:, j]), 0.0, @view(V_ky[:, j]))

            #@timeit to "VecMat" V_kx[j, :] = U_k[j, :]' * simul.DxMatrices[(k-1)*pointsPerElement + j]
            #@timeit to "VecMat" V_kx[j, :] = adjoint(@view U_k[j, :]) * simul.DxMatrices[(k-1)*pointsPerElement + j]
            #@timeit to "VecMat" mul!(@view(V_kx[j, :]), simul.DxMatrices[(k-1)*pointsPerElement + j]', @view(U_k[j, :]))
            #@timeit to "VecMat" mul!(scratch, simul.DxMatrices[(k-1)*pointsPerElement + j], U_k[j, :])
            #@timeit to "VecMat" V_kx[j, :] = scratch

            #@timeit to "MatVec" V_ky[:, j] = simul.DyMatrices[(k-1)*pointsPerElement + j] * U_k[:, j]
            #@timeit to "MatVec" V_ky[:, j] = simul.DyMatrices[(k-1)*pointsPerElement + j] * @view U_k[:, j]
            #@timeit to "MatVec" mul!(@view(V_ky[:, j]), simul.DyMatrices[(k-1)*pointsPerElement + j], @view(U_k[:, j]))
            #@timeit to "MatVec" V_ky[:, j] = U_k[:, j]' * transpose(simul.DyMatrices[(k-1)*pointsPerElement + j])
            #@timeit to "MatVec" mul!(V_ky[:, j], simul.DyMatrices[(k-1)*pointsPerElement + j], @view U_k[:, j])
            #@timeit to "MatVec" mul!(scratch, simul.DyMatrices[(k-1)*pointsPerElement + j], U_k[:, j])
            #@timeit to "MatVec" V_ky[:, j] = scratch

        end
        
    
        @timeit to "weightsMul" V_kx .= simul.QuadWeights .* V_kx
        @timeit to "weightsMul" V_ky .= simul.QuadWeights' .* V_ky
        
        #V_kx = (delta_y_k/delta_x_k) * simul.QuadWeights .* V_kx
        #V_ky = (delta_x_k/delta_y_k) * simul.QuadWeights' .* V_ky

        @timeit to "add" V_k .= V_kx .+ V_ky

        
        @timeit to "set" SetDegreesOfFreedom!(simul, k, laplaceVals, V_k, true)
        #SetDegreesOfFreedom!(simul, k, laplaceVals, V_k_old, true)

        #println("max abs difference between old and new methods: " * string(maximum(abs.(V_k-V_k_old)))) # this difference was something on the order of 1e-15 when I tried it

    end


    @timeit to "adjustments" laplaceVals .= .- (simul.timestep^2 ./ (simul.M .+ 0.5 .* (alpha/beta).*simul.timestep.*simul.M_b)) .* laplaceVals   #u Laplace v = div(u grad v) - grad u \cdot grad v, hence the sign

    #show(to)

    return laplaceVals

end



function ParLaplaceTerm(simul::SEM_Wave, U::Matrix{Float64}) ############# Claude-generated from LaplaceTerm #############

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

    alpha = simul.bc[1]; beta = simul.bc[2]

    if (simul.useMMS == false)
        forcingVals .= forcingVals .* cos(simul.omega*t)
    end

    #println("older version: ")
    #println(-(simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* forcingVals)
    
    return -(simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* forcingVals    # negative sign because we are solving u_tt = u_xx + u_yy - fcos(\omega x) as opposed to u_tt = u_xx + u_yy + fcos(\omega x).

end


function BoundaryTerm(simul::SEM_Wave, stepnumber::Int64)

    boundaryVals = zeros(length(simul.y), length(simul.x)) # not efficient, most entries of the matrix will be zero

    t = simul.timestep*(stepnumber-1) # current timepoint

    alpha = simul.bc[1]; beta = simul.bc[2] # the case where beta = 0 has to be handled separately...
    QuadWeights = simul.QuadWeights
    xNodes = simul.xNodes
    yNodes = simul.yNodes

    pointsPerElement = simul.N + 1                                      # number of quadrature points

    Kx = length(simul.xNodes) - 1                                       # number of elements in the x-direction
    Ky = length(simul.yNodes) - 1                                       # number of elements in the y-direction        


    if simul.useMMS

        # if we use MMS, g should just be whatever the target solution has 
        refNormalDer, refTimeDer = SEM_Wave_2d.BoundaryDersMMS(simul, t)

        simul.g = alpha * refTimeDer + beta .* simul.c_square .* refNormalDer 




    end

    # if no MMS, simul.g already stores the correct information about g.
    # Now integrate g along the boundary
    boundaryVals = simul.M_b .* simul.g

    #println("older version!")
    #println(typeof((1/beta) * (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* boundaryVals))
    #println((1/beta) * (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* boundaryVals )

    return (1/beta) * (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* boundaryVals 

end


function TimeSteppingTerm(simul::SEM_Wave)

    alpha = simul.bc[1]; beta = simul.bc[2]    
    term = 2*simul.M .* simul.uNow + (0.5*(alpha/beta)*simul.timestep*simul.M_b - simul.M) .* simul.uPrev

    #println("older version: ")
    #println(term ./ (simul.M + 0.5*(alpha/beta) * simul.timestep * simul.M_b))

    return term ./ (simul.M + 0.5*(alpha/beta) * simul.timestep * simul.M_b)

end


function TimeSteppingTerm!(simul::SEM_Wave)

    alpha = simul.bc[1]; beta = simul.bc[2]    
    simul.scratch .= 2*simul.M .* simul.uNow .+ (0.5*(alpha/beta)*simul.timestep*simul.M_b - simul.M) .* simul.uPrev

    simul.uNext .= simul.uNext .+ simul.scratch ./ (simul.M + 0.5*(alpha/beta) * simul.timestep * simul.M_b)

end


function Initialise!(simul::SEM_Wave, uStart::Matrix{Float64}, uStartDer::Matrix{Float64}, Tend::Float64, nsteps::Int64, forcing::Matrix{Float64}, omega::Float64, bc::Vector{Float64}, g::Matrix{Float64})

    #time-related things
    simul.Tend = Tend
    simul.nsteps = nsteps
    simul.timestep = Tend/nsteps

    simul.omega = omega
    simul.fVals = forcing

    simul.bc = bc
    simul.g = g

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

    alpha = simul.bc[1]
    beta = simul.bc[2]

    SU = (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b) .* LaplaceTerm(simul, uStart) / simul.timestep^2

    G = zeros(length(simul.y), length(simul.x))


    G = (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b) .* BoundaryTerm(simul, 1) / simul.timestep^2
    F = (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b) .* ForcingTerm(simul, 1) / simul.timestep^2

    simul.uPrev = uStart - simul.timestep*uStartDer + 0.5*simul.timestep^2 * (1 ./ simul.M).*(G - (alpha/beta) * simul.M_b .* uStartDer +
                                                                                            SU +
                                                                                            F)

    
    if (simul.useMMS == true)
        for j = 1:length(simul.x)
            for i = 1:length(simul.y)
                simul.uPrev[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], -simul.timestep, 0, 0, 0, simul.MMS_j)
            end
        end
    end
    

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



function GetDegreesOfFreedom(simul::SEM_Wave, k::Int64, u::Matrix{Float64})
    
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


function SetDegreesOfFreedom!(simul::SEM_Wave, k::Int64, v::Matrix{Float64}, v_k::Matrix{Float64}, add::Bool)

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



function InitialiseMMS(simul::SEM_Wave, uStart::Matrix{Float64}, uStartDer::Matrix{Float64}, Tend::Float64, nsteps::Int64, forcing::Matrix{Float64}, omega::Float64, bc::Vector{Float64}, g::Matrix{Float64})

    Initialise!(simul, uStart, uStartDer, Tend, nsteps, forcing, omega, bc, g)

    uPrevRef = zeros(length(simul.y), length(simul.x))

    for j = 1:length(simul.x)
        for i = 1:length(simul.y)
            uPrevRef[i, j] = MMS.MMSfun(simul.x[j], simul.y[end+1-i], -simul.timestep, 0, 0, 0, simul.MMS_j)
        end
    end

    return [simul.uPrev, uPrevRef]

end


function Waveholtz(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, bc, g, tol::Float64)

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))
    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))
    #nsteps = Integer(ceil(0.5 * Tend * cMax * (1/delta_x + 1/delta_y)))


    timestep = Tend/nsteps

    #println("Number of timesteps for waveholtz: " * string(nsteps) * " with a step size of " * string(timestep))

    # starting guess
    uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = zeros(length(simul.y), length(simul.x))
    oldAppx = zeros(2*length(simul.y), length(simul.x)) # stores uFiltered and uFilteredDer from the previous iteration

    animate = false     # we do not want any animations of the wave eq solutions
    res = Inf


    nIter = 0
    maxIter = 1e5


    while res > tol

        to = TimerOutput()

        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, animate)

        #show(to)

        res = SEM_Wave_2d.SeminormWH(simul, uStart - simul.uFiltered, uStartDer - simul.uDerFiltered) / SEM_Wave_2d.SeminormWH(simul, uStart, uStartDer)

        #res = (SEM_Wave_2d.GradIntegral(simul, uStart - simul.uFiltered, uStart - simul.uFiltered) + SEM_Wave_2d.LpNorm(simul, uStartDer, simul.uDerFiltered, 2)^2)^(1/2) / 
        #        (SEM_Wave_2d.GradIntegral(simul, simul.fVals, simul.fVals))^(1/2)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered
        

        #=
        # relative residual
        if nIter > 1

            res = (LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)^2 + 
                LpNorm(simul, simul.uDerFiltered, oldAppx[length(simul.y) + 1:end, :], 2)^2)^(0.5) / LpNorm(simul, simul.fVals, zeros(length(simul.y), length(simul.y)), 2)

        end
        =#

        oldAppx = [simul.uFiltered; simul.uDerFiltered]
        

        nIter = nIter + 1

        if nIter > maxIter
            res = 0.0 # i.e. break the run
            println("maxIter reached for omega = " * string(omega))
        end

        #println("$omega || iteration: " * string(nIter) * " || residual: " * string(res))

    end

    
    return simul.uFiltered, nIter

end


function Waveholtz(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, bc, g, maxIter::Int64)

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))
    nsteps = Integer(ceil(1.5*Tend * cMax * (1/delta_x + 1/delta_y)))

    # starting guess
    uStart = zeros(length(simul.y), length(simul.x))
    uStartDer = zeros(length(simul.y), length(simul.x))
    oldAppx = zeros(2*length(simul.y), length(simul.x)) # stores uFiltered and uFilteredDer from the previous iteration

    animate = false     # we do not want any animations of the wave eq solutions
    res = NaN
    
    nIter = 0

    while maxIter > nIter
        
        println("iteration: " * string(nIter) * " || residual: " * string(res))

        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, animate)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

        # relative residual
        res = (LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)^2 + 
               LpNorm(simul, simul.uDerFiltered, oldAppx[length(simul.y) + 1:end, :], 2)^2)^(0.5) /
              (LpNorm(simul, oldAppx[1:length(simul.y), :], zeros(length(simul.y), length(simul.x)), 2)^2 + 
               LpNorm(simul, oldAppx[length(simul.y) + 1:end, :], zeros(length(simul.y), length(simul.x)), 2)^2)^0.5
        
        #res = LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)/LpNorm(simul, simul.uFiltered, zeros(length(simul.y), length(simul.y)), 2)

        oldAppx = [simul.uFiltered; simul.uDerFiltered]

        nIter = nIter + 1

    end

    return simul.uFiltered, res

end


function WaveholtzGMRES(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, bc, g, tol::Float64)
    

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))

    nsteps = Integer(ceil(1.1*Tend * cMax * (1/delta_x + 1/delta_y)))

    timestep = Tend/nsteps

    println("time step size for the wave solver: " * string(timestep))
    
    nx = length(simul.x)
    ny = length(simul.y)
    N = Int(length(simul.x) * length(simul.y)) # half of the number of degrees of freedoms of our system
    

    # perfrom one WH step on the zero vecto, this will be the RHS of our linear system
    uStart = zeros(ny, nx)
    uStartDer = zeros(ny, nx)

    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, false)
    b = [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)]
    
    # these should store the iterates, of which there are at most 2N but hopefully much fewer
    Q = spzeros(2*N, 2*N)            #Q = zeros(ComplexF64, n,m+1);
    H = spzeros(2*N, 2*N)            #H = zeros(ComplexF64, m+1,m);
    Q[:,1] = b/norm(b);
    res = Inf

    k = 1

    appx = spzeros(2*N)
    next_appx = spzeros(2*N)
    
    # the matvecs are expensive, so we evaluate the residual only when k mod testFreq = 0.
    resFreq = 20;


    while res > tol
        
        #linear operator L: v \mapsto v - \Pi v + b    (b = \Pi 0)
        
        # reshaping the columns of Q into a matrix and running another wave equation solve
        SEM_Wave_2d.Simulate(simul, Matrix(reshape(Q[1:Int(N), k], ny, nx)), Matrix(reshape(Q[(Int(N)+1):end, k], ny, nx)), Tend, nsteps, fVals, omega, bc, g, false)
        
        w = Q[:, k] - [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)] + b            # Matrix-vector product with last element
        
        # Orthogonalize w against columns of Q.
        h, beta, z = double_gs(Q, w, k);
        
        #Put Gram-Schmidt coefficients into H
        H[1:(k+1), k] = [h; beta];
        
        # normalize
        Q[:, k+1] = z/beta;

        e_1 = zeros(k + 1)
        e_1[1] = 1

        # update variables
        appx = next_appx
        next_appx = Q[:, 1:k]*((H[1:(k+1), 1:k]\e_1)*norm(b))
        k = k+1

        if (k % resFreq == 1)
            SEM_Wave_2d.Simulate(simul, Matrix(reshape(next_appx[1:N], ny, nx)), Matrix(reshape(next_appx[(N+1):end], ny, nx)), Tend, nsteps, fVals, omega, bc, g, false)
            res = norm(next_appx - [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)], 2) / norm(next_appx, 2)
            println("iteration: " * string(k) * " || residual: " * string(res))
        else
            println("iteration: " * string(k) * " || ")
        end

        # some alternatives for the stopping condition that I do not recommend

        #res = norm(appx - next_appx)/norm(appx)
        
        #u0 = Matrix(reshape(appx[1:N, 1], ny, nx))
        #u0_next = Matrix(reshape(next_appx[1:N, 1], ny, nx))
        #res = SEM_Wave_2d.LpNorm(simul, u0_next, u0, 2) / SEM_Wave_2d.LpNorm(simul, u0, zeros(length(simul.y), length(simul.x)), 2)

        #res = ErrorEstimate(simul, u0, fVals, omega)
        

    end
    
    u_0 = Matrix(reshape(next_appx[1:N, 1], ny, nx))
    u_1 = Matrix(reshape(next_appx[N+1:end, 1], ny, nx))

    return u_0, u_1, k

end

function WaveholtzGMRESnew(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, bc, g, tol)

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))

    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))

    timestep = Tend/nsteps

    #println("time step size for the wave solver: " * string(timestep))

    nx = length(simul.x)
    ny = length(simul.y)
    N = Int(length(simul.x) * length(simul.y)) # half of the number of degrees of freedoms of our system




    

    # perfrom one WH step on the zero vector, this will be the RHS of our linear system
    uStart = zeros(ny, nx)
    uStartDer = zeros(ny, nx)

    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, false)
    b = [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)]


    function WaveholtzAction(vec)    

        SEM_Wave_2d.Simulate(simul, Matrix(reshape(vec[1:Int(N)], ny, nx)), Matrix(reshape(vec[(Int(N)+1):end], ny, nx)), Tend, nsteps, fVals, omega, bc, g, false)        
        w = vec - [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)] + b

        return w
        
    end

    WaveholtzMatVec = LinearMap(WaveholtzAction, 2*N) # the matvec as a linear map

    #x, history = gmres(WaveholtzMatVec, b, verbose=true)

    #x, history = gmres(WaveholtzMatVec, b, log=true, verbose=true, reltol=tol,  restart=100)    
    x, history = gmres(WaveholtzMatVec, b, log=true, reltol=tol,  restart=100)    


    u_0 = Matrix(reshape(x[1:N, 1], ny, nx))
    u_1 = Matrix(reshape(x[N+1:end, 1], ny, nx))

    return u_0, u_1, history

end

function WaveholtzConvGMRES(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, bc, g, nIter::Int64)

    ### produces a plot of the convergence history of nIter gmres solves
    

    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))

    nsteps = Integer(ceil(1.1*Tend * cMax * (1/delta_x + 1/delta_y))) 

    timestep = Tend/nsteps

    println("time step size for the wave solver: " * string(timestep))
    
    nx = length(simul.x)
    ny = length(simul.y)
    N = Int(length(simul.x) * length(simul.y)) # half of the number of degrees of freedoms of our system
    

    # perfrom one WH step on the zero vector, this will be the RHS of our linear system
    uStart = zeros(ny, nx)
    uStartDer = zeros(ny, nx)

    SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, false)
    b = [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)]
    
    #these should store the iterates, of which there are at most 2N but hopefully much fewer
    Q = spzeros(2*N, 2*N)            #Q = zeros(ComplexF64, n,m+1);
    H = spzeros(2*N, 2*N)            #H = zeros(ComplexF64, m+1,m);
    Q[:,1] = b/norm(b);
    res = Inf

    k = 1

    appx = spzeros(2*N)
    next_appx = spzeros(2*N)

    convHist = zeros(nIter, 1)


    while nIter > k
        
        #linear operator L: v \mapsto v - \Pi v + b (b = \Pi 0)
        #SEM_Wave_1d.Simulate(simul, Vector(Q[1:Int(dofs/2), k]), Vector(Q[(Int(dofs/2)+1):end, k]), false)
        
        # reshaping the column of Q into a matrix and running another wave equation solve
        SEM_Wave_2d.Simulate(simul, Matrix(reshape(Q[1:Int(N), k], ny, nx)), Matrix(reshape(Q[(Int(N)+1):end, k], ny, nx)), Tend, nsteps, fVals, omega, bc, g, false)
        
        w = Q[:, k] - [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)] + b            # Matrix-vector product with last element
        
        # Orthogonalize w against columns of Q.
        h, beta, z = double_gs(Q, w, k);
        
        #Put Gram-Schmidt coefficients into H
        H[1:(k+1), k] = [h; beta];
        
        # normalize
        Q[:, k+1] = z/beta;

        e_1 = zeros(k + 1)
        e_1[1] = 1

        #update variables
        appx = next_appx
        next_appx = Q[:, 1:k]*((H[1:(k+1), 1:k]\e_1)*norm(b))
        k = k+1

        # res = L(next_appx) - b = (I - \Pi)next_appx
        SEM_Wave_2d.Simulate(simul, Matrix(reshape(next_appx[1:N], ny, nx)), Matrix(reshape(next_appx[(N+1):end], ny, nx)), Tend, nsteps, fVals, omega, bc, g, false)
        res = norm(next_appx - [reshape(simul.uFiltered, N, 1); reshape(simul.uDerFiltered, N, 1)], 2) / norm(next_appx, 2)
        #res = norm(appx - next_appx)/norm(appx)
        #=
        u0 = Matrix(reshape(appx[1:N, 1], ny, nx))
        u0_next = Matrix(reshape(next_appx[1:N, 1], ny, nx))
        res = SEM_Wave_2d.LpNorm(simul, u0_next, u0, 2) / SEM_Wave_2d.LpNorm(simul, u0, zeros(length(simul.y), length(simul.x)), 2)
        =#
        #res = ErrorEstimate(simul, u0, fVals, omega)
        println("iteration: " * string(k) * " || residual: " * string(res))
        
        convHist[k] = res

    end
    
    u_0 = Matrix(reshape(next_appx[1:N, 1], ny, nx))
    u_1 = Matrix(reshape(next_appx[N+1:end, 1], ny, nx))
    #return next_appx, LpNorm(simul, Vector(next_appx), zeros(length(simul.x)), 2)

    plt = plot(2:nIter, convHist[2:end], yscale=:log10)
    savefig(plt, "gmresConvHistory")

    return u_0, u_1, convHist

end

function WaveholtzAnimation(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, bc, g, maxIter::Int64, logscale = false)

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
        
        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, animate)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered
        #uStartDer = zeros(length(simul.y), length(simul.x))


        # relative residual
        res = (LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.y), :], 2)^2 + 
               LpNorm(simul, simul.uDerFiltered, oldAppx[length(simul.y) + 1:end, :], 2)^2)^(0.5) /
              (LpNorm(simul, oldAppx[1:length(simul.y), :], zeros(length(simul.y), length(simul.x)), 2)^2 + 
               LpNorm(simul, oldAppx[length(simul.y) + 1:end, :], zeros(length(simul.y), length(simul.x)), 2)^2)^0.5
        

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

function WaveholtzConvHistory(simul::SEM_Wave, omega::Float64, fVals::Matrix{Float64}, bc, g, maxIter::Int64, makePlot::Bool = true)

    ### produces a plot of the convergence history of nIter Waveholtz iterations

    # plots the convergence history of a Waveholtz iteration
    # Waveholtz-appropriate parameters for the wave solver
    Tend = 2*pi/omega
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    delta_y = minimum(simul.y[2:end] - simul.y[1:end-1])
    cMax = sqrt(maximum(simul.c_square))
    #nsteps = Integer(ceil(1.5*Tend * cMax * (1/delta_x + 1/delta_y)))
    nsteps = Integer(ceil(Tend * cMax * (1/delta_x + 1/delta_y)))

    nsteps = nsteps*2
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

        SEM_Wave_2d.Simulate(simul, uStart, uStartDer, Tend, nsteps, fVals, omega, bc, g, animate)

        #fValsH1Norm = (SEM_Wave_2d.LpNorm(simul, simul.fVals, zeros(length(simul.y), length(simul.x)), 2)^2 + SEM_Wave_2d.GradIntegral(simul, simul.fVals, simul.fVals))^(1/2)
        #res = (SEM_Wave_2d.GradIntegral(simul, uStart - simul.uFiltered, uStart - simul.uFiltered) + SEM_Wave_2d.LpNorm(simul, uStartDer, simul.uDerFiltered, 2)^2)^(1/2) / fValsH1Norm


        res = SEM_Wave_2d.SeminormWH(simul, uStart - simul.uFiltered, uStartDer - simul.uDerFiltered) / SEM_Wave_2d.SeminormWH(simul, uStart, uStartDer)

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

function ErrorEstimate(simul::SEM_Wave, u::Matrix{Float64}, fVals::Matrix{Float64}, omega::Float64) 
    
    ### should estimate || Laplace u + omega^2 u - f ||_{interior nodes} + || alpha u + \beta \omega \partial_n u ||_{boundary}, which should show whether u is a good approximation to the corresponding PDE solution

    println("ErrorEstimate is under construction: do not use it unless you understand the risks")

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

        integrand_k .= (abs.(SEM_Wave_2d.GetDegreesOfFreedom(simul, k, u) .- SEM_Wave_2d.GetDegreesOfFreedom(simul, k, v))).^p

        
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
    out = (abs(SEM_Wave_2d.GradIntegral(simul, U_real, U_real) + SEM_Wave_2d.LpNorm(simul, U_imag, zeros(length(simul.y), length(simul.x)), 2)^2))^(0.5)

    return out

end


function GradIntegral(simul::SEM_Wave, U::Matrix{Float64}, V::Matrix{Float64})
    ### calculates \int_\Omega c^2(x, y) \grad u \cdot \grad v dxdy where the nodal values of u and v are stored in the matrices U and V
    ### a generalisation of LaplaceTerm


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
        c_square_k = GetDegreesOfFreedom(simul, k, simul.c_square)

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
    if abs(gradIntegral) < 1e-16 # importantly, if there were a bug that made gradIntegral < -0.1, say, this would not hide that problem.
        gradIntegral = 0.0
    end

    return gradIntegral

end

function SobolevNorm(simul::SEM_Wave, U::Matrix{Float64})

    ### calculates the H1-norm of u whose values is stored in the matrix U
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

    return LpNorm(simul, U, zeros(length(simul.y), length(simul.x)), 2) + sqrt(gradIntegral)


end


function BoundaryDersMMS(simul::SEM_Wave, t::Float64)
    ### calculates the normal and time derivatives of the MMS solution on the boundary at time t for use in MMS tests

    refNormalDer = zeros(length(simul.y), length(simul.x))
    refTimeDer = zeros(length(simul.y), length(simul.x))

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


    # we need to take a weighted average in order to take account for the cusp at the corners
    delta_x = simul.xNodes[2] - simul.xNodes[1]
    delta_y = simul.yNodes[2] - simul.yNodes[1]

        
    
    refNormalDer[1, 1] = (-MMS.MMSfun(simul.x[1], simul.y[end], t, 1, 0, 0, simul.MMS_j)*delta_y + MMS.MMSfun(simul.x[1], simul.y[end], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 
    refNormalDer[1, end] = (MMS.MMSfun(simul.x[end], simul.y[end], t, 1, 0, 0, simul.MMS_j)*delta_y + MMS.MMSfun(simul.x[end], simul.y[end], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 
    refNormalDer[end, 1] = (-MMS.MMSfun(simul.x[1], simul.y[1], t, 1, 0, 0, simul.MMS_j)*delta_y - MMS.MMSfun(simul.x[1], simul.y[1], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 
    refNormalDer[end, end] = (MMS.MMSfun(simul.x[end], simul.y[1], t, 1, 0, 0, simul.MMS_j)*delta_y - MMS.MMSfun(simul.x[end], simul.y[1], t, 0, 1, 0, simul.MMS_j)*delta_x)/(delta_x + delta_y) 
    

    return refNormalDer, refTimeDer

end



end