
### mms.jl has changed since this was written, so if the mms functionality is going to be used, changes will be needed. cf. SEM_Wave_2d.jl ###

module SEM_Wave_1d

using LinearAlgebra
using SparseArrays
using Plots
using FastGaussQuadrature
using Polynomials
using TimerOutputs

include("MMS.jl")
using .MMS


export SEM_Wave

mutable struct SEM_Wave
    #stores data for SEM approximation of the wave equation \frac{\partial^2 u}{\partial t^2} = \frac{\partial^2 u}{\partial x^2} - f(x)\cos(\omega t)


    nodes::Vector{Float64}          # the endpoints of the elements
    N::Int64                        # the degree of the Lagrange polynomials; N = pointsPerElement + 1
    x::Vector{Float64}              # the actual coordinates between nodes[1] and nodes[end] where we approximate our solution
    c_square::Function              # function that gives the squared wave speed in a point

    nsteps::Int64                   # number of steps we make in time
    Tend::Float64                   # the endpoint of our interval in time
    timestep::Float64               # size of the timesteps (= Tend/nsteps)
    
    uPrev::Vector{Float64}          # value of the displacement in the previous timepoint
    uNow::Vector{Float64}           # value of the displacement in the current timepoint
    uNext::Vector{Float64}          # value of the displacement in the next timepoint
    uFiltered::Vector{Float64}      # the filtered time average of u for use in Waveholtz
    uDerFiltered::Vector{Float64}   # the filtered time average of u_t for use in Waveholtz
    useWaveholtz::Bool              # decides whether to filter the solution of the wave equation or not

    fVals::Vector{Float64}          # values of the driving term of the equation
    omega::Float64                  # frequency parameter omega
    bc::Vector{Float64}             # contains the \alpha and \beta used in the bc 
    g::Vector{Float64}              # function values for g in Boundary condition
    
    
    QuadPoints::Vector{Float64}     # the points in [-1, 1] we use in our Gauss-Lobatto quadrature
    QuadWeights::Vector{Float64}    # the weights corresponding to the points in QuadPoints
    G::Array{Float64}               # the reference stiffness matrix: G_{jm} = \sum_{l=0}^{N} w_l l'_m(\xi_l) l'_j(\xi_l); w_l \in QuadWeights, \xi_l \in QuadPoints, l_j's are the Lagrange polynomials                            
    M::Vector{Float64}              # the mass matrix, which is diagonal => stored as vector
    M_b::Vector{Float64}            # the mass matrix on the boundary, used to handle terms from the boundary intergral
    
    useMMS::Bool                    # determines whether an MMS test is run or not
    MMS_j::MMS_jet                  # stores the relevant data for the MMS test


    function SEM_Wave(
        boundaryCoordinates::Vector{Float64},   # coordinates of xl and xr
        K::Int64,                               # number of elements
        N::Int64,                               # degree of interpolation polynomials: N + 1 points in every element
        c_square::Function                      # wave speed
        )
        
        QuadPoints, QuadWeights = gausslobatto(N + 1)
        uPrev = zeros(K*N + 1)
        uNow = zeros(K*N + 1)
        uNext = zeros(K*N + 1)
        uFiltered = zeros(K*N + 1)
        uDerFiltered = zeros(K*N + 1)

        nodes = collect(LinRange(boundaryCoordinates[1], boundaryCoordinates[2], K+1))
        x = ConstructX(nodes, QuadPoints)

        fVals = zeros(K*N + 1)
        omega = 0.0
        bc = [0.0, 1.0]
        g = [0.0, 0.0]
        
        
        useMMS = false
        useWaveholtz = false

        Tend = 1                # arbitrary values, just here as placeholders
        nsteps = 100
        timestep = Tend/nsteps

        # make G and the inverse of the mass matrix
        G = ConstructG(QuadPoints, QuadWeights)
        M, M_b = ConstructMs(nodes, QuadWeights)

        # change some coefficients in the MMS
        MMS_j = MMS.MMS_jet(ones(2, 1), 2)
        MMS_j.coeff[1, 2] = 2.0
        MMS_j.coeff[2, 1] = exp(1)
        MMS_j.coeff[2, 3] = -2
        MMS_j.coeff[1, 3] = (1 + sqrt(5))/2
        
                new(nodes, N, x, c_square, nsteps, Tend, timestep,
                uPrev, uNow, uNext, uFiltered, uDerFiltered, useWaveholtz, fVals, omega, bc, g, QuadPoints, QuadWeights, G, M, M_b, useMMS, MMS_j)
        
    end

end


function Simulate(simul::SEM_Wave, uStart, uStartDer, Tend::Float64, nsteps::Int64, forcing::Vector{Float64}, omega::Float64, bc::Vector{Float64}, g::Vector{Float64}, animate = false, snapshotFrequency = 10, plotHeight = 10.0, animationName = "test.gif")
    # approximates a solution to the wave equation using SEM_Wave
    # creates a gif named animationName if animate == true
    # snapshotFrequency decides how many steps are made before a frame is saved in the gif
    # the y-axis shown in the gif is [-plotHeight, plotheight]

    plotheight = maximum(uStart)

    Initialise!(simul, uStart, uStartDer, Tend, nsteps, forcing, omega, bc, g)
    
    #create an Animation if one is asked for
    if (animate == true)
       
        anim = Animation()

        plot(simul.x, simul.uNow, xlims=(simul.nodes[1], simul.nodes[end]), ylims=(-plotHeight, plotHeight), legend=:false)
        frame(anim)    
    
    end
    
    
    #make the steps
    for n = 1:simul.nsteps


        MakeStep!(simul, n)    

        #making the Animation is very slow compared to the actual calculations made here.
        #store a frame in our Animation once every snapshotFrequency steps
        if (animate == true && n%snapshotFrequency == 0)
            #plot(simul.x, simul.uNow, xlims=(nodes[1],nodes[end]), ylims=(-0.02, 0.02))
            plot(simul.x, simul.uNow, xlims=(simul.nodes[1], simul.nodes[end]), ylims=(-plotHeight, plotHeight), legend=:false)
            frame(anim)
            #println(n)
        end

            
    end

    #finally make and store the gif
    if (animate == true)

        gif(anim, animationName, fps=40) #if you need the gif to run a different framerate, change the value of fps

    end

    
    #return simul.uNow
    
end


function Waveholtz(simul::SEM_Wave, forcing::Vector{Float64}, omega::Float64, bc::Vector{Float64}, g::Vector{Float64}, tol::Float64)

    # Waveholtz-specific parameters for the wave solver
    simul.fVals = -forcing
    simul.omega = omega
    simul.Tend = 2*pi/omega
    simul.bc = bc

    # let simul know that we are using Waveholtz and want the filtered solution
    simul.useWaveholtz = true

    # appropriate parameters for the wave solver:
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    nsteps = Integer(ceil(1.5*simul.Tend/delta_x))


    # starting guess
    uStart = zeros(length(simul.x))
    uStartDer = zeros(length(simul.x))
    oldAppx = zeros(2*length(simul.x))

    animate = false     #we do not want any animations of the wave eq solutions
    res = Inf


    nIter = 0
    maxIter = 20000

    while res > tol

        SEM_Wave_1d.Simulate(simul::SEM_Wave, uStart, uStartDer, simul.Tend::Float64, nsteps::Int64, simul.fVals::Vector{Float64}, omega::Float64, bc::Vector{Float64}, g::Vector{Float64})

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

        #println(maximum(uStart))

        #println(norm(simul.uFiltered - oldAppx[1:length(simul.x)]))

        #relative H2-residual
        res = (LpNorm(simul, simul.uFiltered, oldAppx[1:length(simul.x)], 2)^2 + 
               LpNorm(simul, simul.uDerFiltered, oldAppx[length(simul.x) + 1:end], 2)^2)^(0.5)/
              (LpNorm(simul, oldAppx[1:length(simul.x)], zeros(length(simul.x)), 2)^2 + 
               LpNorm(simul, oldAppx[length(simul.x) + 1:end], zeros(length(simul.x)), 2)^2)^0.5

        oldAppx = [simul.uFiltered; simul.uDerFiltered]
        
        println(res)
        
        nIter = nIter + 1

        #println(nIter)

        if nIter > maxIter
            res = 0.0 # i.e. break the run
            println("maxIter reached for omega = " * string(omega))
        end

    end

    #return simul.uFiltered, LpNorm(simul, simul.uFiltered, zeros(length(simul.x)), 2)

    #no longer interested in filtered solutions:
    simul.useWaveholtz = false

    return simul.uFiltered, nIter

end


function WaveholtzGMRES(simul::SEM_Wave, forcing::Vector{Float64}, omega::Float64, bc::Vector{Float64}, g::Vector{Float64}, tol::Float64)
    

    #Waveholtz-specific parameters for the wave solver
    simul.fVals = -forcing
    simul.omega = omega
    simul.Tend = 2*pi/omega
    simul.bc = bc
    
    #let simul know that we are using waveholtz and want the filtered solution
    simul.useWaveholtz = true
    animate = false     #we do not want any animations of the wave eq solutions

    #appropriate parameters for the wave solver:
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    nsteps = Integer(ceil(1.5*simul.Tend/delta_x))


    # perform one WH step on the zero vector 
    zeroVec = zeros(length(simul.x))

    SEM_Wave_1d.Simulate(simul::SEM_Wave, zeroVec, zeroVec, simul.Tend::Float64, nsteps::Int64, simul.fVals::Vector{Float64}, omega::Float64, bc::Vector{Float64}, g::Vector{Float64})
    b = [simul.uFiltered; simul.uDerFiltered]
    n = length(b);
    
    #these should store the iterates, of which there are at most n but hopefully much fewer
    Q = spzeros(n, n)            #Q = zeros(ComplexF64, n,m+1);
    H = spzeros(n, n)            #H = zeros(ComplexF64, m+1,m);
    Q[:,1] = b/norm(b);
    res = Inf

    k = 1

    appx = spzeros(n)
    next_appx = spzeros(n)


    while res > tol
        
        #linear operator L: v \mapsto v - \Pi v + b (b = \Pi 0)
        SEM_Wave_1d.Simulate(simul, Vector(Q[1:Int(n/2), k]), Vector(Q[(Int(n/2)+1):end, k]), simul.Tend::Float64, nsteps::Int64, simul.fVals::Vector{Float64}, omega::Float64, bc::Vector{Float64}, g::Vector{Float64})
        w = Q[:, k] - [simul.uFiltered; simul.uDerFiltered] + b            # Matrix-vector product with last element
        
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

        res = LpNorm(simul, Vector(appx), Vector(next_appx), 2)/LpNorm(simul, Vector(appx), zeros(length(appx)), 2)
        println(res)

    end

    simul.useWaveholtz = false
    
    #return next_appx, LpNorm(simul, Vector(next_appx), zeros(length(simul.x)), 2)
    return next_appx, k

end


#=
function Waveholtz(simul::SEM_Wave, omega::Float64, fVals::Vector{Float64}, nIter::Int64) 

    #Waveholtz-specific parameters for the wave solver
    simul.fVals = -fVals
    simul.omega = omega
    simul.Tend = 2*pi/omega
    

    #appropriate parameters for the wave solver:
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
    simul.timestep = simul.Tend/simul.nsteps

    #starting guess
    uStart = zeros(length(simul.x))
    uStartDer = zeros(length(simul.x))

    
    animate = false     #we do not want any animations of the wave eq solutions

    for j = 1:nIter
        
        SEM_Wave_1d.Simulate(simul, uStart, uStartDer, animate)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

    end

end
=#

function Initialise!(simul::SEM_Wave, uStart::Vector{Float64}, uStartDer::Vector{Float64}, Tend::Float64, nsteps::Int64, forcing::Vector{Float64}, omega::Float64, bc::Vector{Float64}, g::Vector{Float64})
    
    # reset Waveholtz
    if simul.useWaveholtz
        simul.uFiltered = zeros(length(simul.x))
        simul.uDerFiltered = zeros(length(simul.x))
    end

    # set up time-related things
    simul.Tend = Tend
    simul.nsteps = nsteps
    simul.timestep = Tend/nsteps


    # set up the forcing and bonudnary conditions, perhaps could be done more nicely when one uses waveholtz...
    simul.fVals = forcing
    simul.omega = omega
    simul.bc = bc
    simul.g = g


    # set up initial and boundary conditions
    simul.uNow = uStart
    alpha = simul.bc[1]
    beta = simul.bc[2]

    # approximate simul.uPrev
    SU = (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b).*LaplaceTerm(simul, uStart) / simul.timestep^2
    
    G = zeros(length(simul.x))

    if simul.useMMS == false
        G[1] = simul.g[1]
        G[end] = simul.g[end]
    else
        G[1] = alpha*MMS.MMSfun(simul.nodes[1], 0.0, 0, 1, simul.MMS_j) -
               beta*MMS.MMSfun(simul.nodes[1], 0.0, 1, 0, simul.MMS_j)

        G[end] = alpha*MMS.MMSfun(simul.nodes[end], 0.0, 0, 1, simul.MMS_j) +
                 beta*MMS.MMSfun(simul.nodes[end], 0.0, 1, 0, simul.MMS_j)
    end


    F = (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b).*ForcingTerm(simul, 1) / simul.timestep^2    

    simul.uPrev = uStart - simul.timestep*uStartDer + 0.5*simul.timestep^2*(1 ./ simul.M).*(simul.M_b.*((1/beta)*G - (alpha/beta)*uStartDer) + 
                                                                                            SU + 
                                                                                            F)

end




function MakeStep!(simul::SEM_Wave, stepnumber::Int64)    

    t = stepnumber*simul.timestep # the time at uNow after(!) this step

    simul.uNext = TimeSteppingTerm(simul) + 
                  LaplaceTerm(simul, simul.uNow) + 
                  ForcingTerm(simul, stepnumber) + 
                  BoundaryTerm(simul, stepnumber)
    
    # use a second-order approximation of the derivative at t using known values
    uDer = ((3/2)*simul.uNext - 2*simul.uNow + 0.5*simul.uPrev) ./ simul.timestep 
    #uDer = (simul.uNext - simul.uNow)./simul.timestep
    uDerPrev = (simul.uNext - simul.uPrev) ./ (2*simul.timestep)
    #uDerPrev = (simul.uNow - simul.uPrev) ./ simul.timestep
    

    ###########
    #update the variables
    simul.uPrev = simul.uNow
    simul.uNow = simul.uNext    


    if (simul.useWaveholtz)

        #filter the solution for use in Waveholtz:
        Know = (cos(simul.omega*t) - 0.25) * (2/simul.Tend)
        Kprev = (cos(simul.omega * (t - simul.timestep)) - 0.25) * (2/simul.Tend)
        
        simul.uFiltered = simul.uFiltered + (simul.uPrev*Kprev + simul.uNow*Know) .* (simul.timestep/2)
        simul.uDerFiltered = simul.uDerFiltered + (uDerPrev*Kprev + uDer*Know) .* (simul.timestep/2)
        #note how this makes the 1 2 ... 2 1 pattern from the trapezoid method appear

    end

end


function LaplaceTerm(simul::SEM_Wave, u::Vector{Float64})

    #to = TimerOutput()
    K = length(simul.nodes)-1      #number of elements
    degreesOfFreedom = K*simul.N + 1
   
    laplaceVals = zeros(degreesOfFreedom)

    u_k = zeros(simul.N + 1)
    v_k = zeros(simul.N + 1)
    G = simul.G

    for k = 1:K

        #=
        @timeit to "delta_x_k" delta_x_k = simul.nodes[k+1] - simul.nodes[k]

        @timeit to "get" u_k .= GetDegreesOfFreedom(simul, k, u)

        #@timeit to "v_k" v_k .= (2/delta_x_k) * (simul.G*u_k)

        @timeit to "v_k" mul!(v_k, simul.G, u_k, 2/delta_x_k, 0)

        @timeit to "Set" SetDegreesOfFreedom!(simul, k, laplaceVals, v_k, true)
        =#

        delta_x_k = simul.nodes[k+1] - simul.nodes[k]

        u_k .= GetDegreesOfFreedom(simul, k, u)
        x_k = simul.x[((k-1)*simul.N + 1):(k*simul.N + 1)]

        #what exactly does this multiply and add again? v_k = (2/delta_x_k) * G * u_k + 0 perhaps? 
        #mul!(v_k, G, u_k, 2/delta_x_k, 0)

        
        mul!(v_k, G, u_k.*simul.c_square.(x_k), 2/delta_x_k, 0)

        SetDegreesOfFreedom!(simul, k, laplaceVals, v_k, true)

    end

    #show(to)    
    alpha = simul.bc[1]; beta = simul.bc[2]

    laplaceVals = - (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* laplaceVals   #u Laplace v = div(u grad v) - grad u \cdot grad v, hence the sign

    return laplaceVals

end


function ForcingTerm(simul::SEM_Wave, stepnumber::Int64)

    N = simul.N                     # degree of interpolation polynomials
    K = length(simul.nodes)-1       # number of elements
    degreesOfFreedom = K*N + 1
   
    forcingVals = zeros(degreesOfFreedom)

    F_k = zeros(N + 1)
    vals_k = zeros(N + 1)
    
    if simul.useMMS

        t = simul.timestep*(stepnumber-1)

        for j = 1:length(simul.x)
            simul.fVals[j] = MMS.MMSfun(simul.x[j], t, 0, 2, simul.MMS_j) - MMS.MMSfun(simul.x[j], t, 2, 0, simul.MMS_j)
        end

    end


    for k = 1:K

        delta_x_k = simul.nodes[k+1] - simul.nodes[k]

        F_k .= GetDegreesOfFreedom(simul, k, simul.fVals)

        vals_k = 0.5*delta_x_k * (simul.QuadWeights.*F_k)

        SetDegreesOfFreedom!(simul, k, forcingVals, vals_k, true)

    end


    alpha = simul.bc[1]; beta = simul.bc[2]


    if (simul.useMMS == true)
        forcingVals = (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* forcingVals
    else
        forcingVals = (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* forcingVals * cos(simul.omega*(stepnumber-1)*simul.timestep)
    end

    return forcingVals

end


function BoundaryTerm(simul::SEM_Wave, stepnumber::Int64)

    boundaryVals = zeros(length(simul.x))

    t = simul.timestep*(stepnumber-1) #current timepoint

    alpha = simul.bc[1]; beta = simul.bc[2]

    if simul.useMMS == false

        boundaryVals[1] = simul.g[1]*simul.c_square(simul.x[1])
        boundaryVals[end] = simul.g[end]*simul.c_square(simul.x[end])

    else    #if we use MMS, the gs should just be whatever the target solution has 

        boundaryVals[1] = alpha*MMS.MMSfun(simul.nodes[1], t, 0, 1, simul.MMS_j) - 
                          beta*MMS.MMSfun(simul.nodes[1], t, 1, 0, simul.MMS_j)

        boundaryVals[end] = alpha*MMS.MMSfun(simul.nodes[end], t, 0, 1, simul.MMS_j) +
                            beta*MMS.MMSfun(simul.nodes[end], t, 1, 0, simul.MMS_j)

    end

    return (1/beta) * (simul.timestep^2 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* boundaryVals 

end


function TimeSteppingTerm(simul::SEM_Wave)
    
    alpha = simul.bc[1]; beta = simul.bc[2]    
    term = 2*simul.M .* simul.uNow + (0.5*(alpha/beta)*simul.timestep*simul.M_b - simul.M) .* simul.uPrev

    return (1 ./ (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b)) .* term

end


function ConstructX(nodes::Vector{Float64}, QuadPoints::Vector{Float64})
    
    K = length(nodes)-1
    N = length(QuadPoints)

    x = zeros(K*(N-1) + 1)

    #puts the correct values in x
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

    #this is where values of c^2(x_i) are most easily added
    return transpose(D)*(QuadWeights.*D)

end


function BarycentricWeights(points::Vector{Float64})

    #returns the barycentric weights of a vector of points; used in calculating the derivative matrix of the Lagrange polynomials

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


function ConstructMs(nodes::Vector{Float64}, QuadWeights::Vector{Float64})

    #constructs the mass matrix and the boundary mass matrix

    K = length(nodes)-1
    N = length(QuadWeights)

    M = zeros(K*(N-1) + 1)
    M_b = zeros(K*(N-1) + 1)
    
    for k = 1:K
            
        start = (k-1)*(N-1)

        delta_x_k = nodes[k+1] - nodes[k]

        M[(start + 1):(start + N)] += 0.5*QuadWeights*delta_x_k
        
    end

    M_b[1] = 1; M_b[end] = 1

    return M, M_b

end


function GetDegreesOfFreedom(simul::SEM_Wave, k::Int64, u::Vector{Float64})
    # returns all of the of entries of u corresponding to element k
    start = (k-1)*(simul.N)
    return u[(start + 1):(start + simul.N + 1)]

end


function SetDegreesOfFreedom!(simul::SEM_Wave, k::Int64, v::Vector{Float64}, v_k::Vector{Float64}, add::Bool)
    # modifies all of the of entries of u corresponding to element k

    start = (k-1)*(simul.N)

    if add
        #v[(start + 1):(start + (simul.N + 1))] .+= v_k
        view(v, (start + 1):(start + (simul.N + 1))) .+= v_k
    else
        v[(start + 1):(start + (simul.N + 1))] .= v_k
    end

end


function LaplaceMMS(simul::SEM_Wave)   

    n = length(simul.x)
    ref = zeros(n)
    refLap = zeros(n) 

    for j = 1:n
        ref[j] = MMS.MMSfun(simul.x[j], 0, 1, simul.MMS_j)
        refLap[j] = MMS.MMSfun(simul.x[j], 2, 1, simul.MMS_j)
    end

    alpha = simul.bc[1]; beta = simul.bc[2]
    #laplaceVals = (simul.M + 0.5*(alpha/beta)*simul.timestep*simul.M_b) .* SEM_Wave_1d.LaplaceTerm(simul, ref) / simul.timestep^2 # why not this?
    laplaceVals = SEM_Wave_1d.LaplaceTerm(simul, ref) / simul.timestep^2

    return [laplaceVals[2:end-1], refLap[2:end-1]] #compare everywhere except on the boundary
    
end


function ForcingMMS(simul::SEM_Wave, stepnumber::Int64)   

    # only works if (simul.useMMS == true)

    t = (stepnumber-1)*simul.timestep
    n = length(simul.x)
    fRef = zeros(n) 

    for j = 1:n
        fRef[j] = MMS.MMSfun(simul.x[j], t, 0, 2, simul.MMS_j) - MMS.MMSfun(simul.x[j], t, 2, 0, simul.MMS_j)
    end

    alpha = simul.bc[1]; beta = simul.bc[2]
    fVals = SEM_Wave_1d.ForcingTerm(simul, stepnumber) / simul.timestep^2

    return [fVals, fRef]
    
end


function InitialiseMMS(simul::SEM_Wave, uStart::Vector{Float64}, uStartDer::Vector{Float64})

    Initialise!(simul, uStart, uStartDer)
    
    uPrevRef = zeros(length(simul.x))
    
    for j = 1:length(simul.x)
        uPrevRef[j] = MMS.MMSfun(simul.x[j], -simul.timestep, 0, 0, simul.MMS_j)
    end

    return [simul.uPrev, uPrevRef]

end


function SetElements(newNodes::SEM_Wave, uStart::Vector{Float64})

    #this for now only works if length(newNodes) == length(simul.nodes)
    if (length(newNodes) == length(simul.nodes))
        simul.nodes = newNodes
        x = ConstructX(simul.nodes, simul.QuadPoints)
        M, M_b = ConstructMs(nodes, QuadWeights)

    else
        println("the new nodes would change the number of elements")
    end

end

################ under construction ################
function HelmholtzMatrix(simul::SEM_Wave, omega::Float64, bc, fVals)
    # returns the matrix corresponding to the discrete Helmholtz operator

    # setup boundary conditions and frequency
    simul.bc = bc
    simul.omega = omega

    nx = length(simul.x)

    Kx = length(simul.nodes) - 1

    H = spzeros(nx, nx)

    e_j = zeros(nx, 1)
    Le_j = zeros(nx, 1)        


    for j = 1:nx

        prevIndex = maximum([j-1, 1])
        e_j[prevIndex] = 0.0
        e_j[j] = 1.0

        kList = []

        if (mod(j, simul.N) == 1) # if on the edge of an element
            
            if (j == 1) # if the first element
                
                kList = [kList; 1]

            elseif (j == nx) # if the final element
                
                kList = [kList; Kx]

            else # else on the edge between two elements

                index1 = Int((j-1 - mod(j-1, simul.N)) / simul.N) + 1
                index2 = index1 - 1

                kList = [kList; index1; index2]
            
            end

        end
        
        # reset matrix for storing the data 
        Le_j .= zeros(size(e_j))        
        
        

        for s = 1:length(kList)

            k = kList[s]
            SEM_Wave_1d.LocalLaplace!(simul, k, e_j, Le_j) # compute action of Laplacian, store the result in Le_j

        end

        # store data
        H[:, j] .= Le_j

        println("generated row $j out of $nx")

    end

    M = spdiagm(simul.M)
    B = spdiagm((simul.bc[1]/simul.bc[2]) * simul.M_b)
    H = (simul.omega^2 * M + (1im*simul.omega)*B + H)       # construct the Helmholtz operator
    
    F = M*fVals;                                            # the matrix F

    return H, F

end


function LocalLaplace!(simul::SEM_Wave, k::Int64, U::Matrix{Float64}, laplaceVals::Vector{Float64})
    # approximates the Laplacian of U in element k and stores the result in laplaceVals
    # intended to replace the corresponding part of the code in LaplaceTerm by a function call

    u_k = zeros(simul.N + 1)
    v_k = zeros(simul.N + 1)
    
    G = simul.G
    delta_x_k = simul.nodes[k+1] - simul.nodes[k]
    u_k .= GetDegreesOfFreedom(simul, k, U)
    x_k = simul.x[((k-1)*simul.N + 1):(k*simul.N + 1)]
    
    mul!(v_k, G, u_k.*simul.c_square.(x_k), 2/delta_x_k, 0)

    SetDegreesOfFreedom!(simul, k, laplaceVals, v_k, true)


end


function WaveholtzAnimation(simul::SEM_Wave, omega::Float64, fVals::Vector{Float64}, nIter::Int64)

    #Waveholtz-specific parameters for the wave solver
    simul.fVals = -fVals
    simul.omega = omega
    simul.Tend = 2*pi/omega
    

    #appropriate parameters for the wave solver:
    delta_x = minimum(simul.x[2:end] - simul.x[1:end-1])
    simul.nsteps = Integer(ceil(1.5*simul.Tend/delta_x))
    simul.timestep = simul.Tend/simul.nsteps

    #starting guess
    uStart = zeros(length(simul.x))
    uStartDer = zeros(length(simul.x))
    
    
    animate = false     #we do not want any animations of the wave eq solutions
    
    anim = Animation()

    plot(simul.x, simul.uFiltered, xlims=(simul.nodes[1], simul.nodes[end]), ylims=(-1, 1), legend=:false)
    for j = 1:nIter

        println(j)
        
        SEM_Wave_1d.Simulate(simul, uStart, uStartDer, animate)

        uStart = simul.uFiltered
        uStartDer = simul.uDerFiltered

        
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
        

        plot(simul.x, simul.uFiltered, xlims=(simul.nodes[1], simul.nodes[end]), ylims=(plotLow, plotHigh), legend=:false)
        frame(anim) 

    end

    gif(anim, "WaveholtzGif.gif", fps=40)

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


function LpNorm(simul::SEM_Wave, u::Vector{Float64}, v::Vector{Float64}, p::Int64)

    weights = simul.QuadWeights

    K = length(simul.nodes)-1 #number of elements
    norm = 0.0

    for k = 1:K

        k = Int(k)

        u_k = SEM_Wave_1d.GetDegreesOfFreedom(simul, k, u)
        v_k = SEM_Wave_1d.GetDegreesOfFreedom(simul, k, v)

        xl = simul.nodes[k]

        xr = simul.nodes[k+1]
        
        #integrate
        norm = norm + 0.5*(xr-xl)*dot(weights, abs.((u_k - v_k).^p))

    end

    norm = norm^(1/p)

    return norm

end


end

