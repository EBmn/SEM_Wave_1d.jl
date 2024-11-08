

module SEM_Wave_1d

using LinearAlgebra
using Plots
using FastGaussQuadrature
using Polynomials
using TimerOutputs

include("MMS.jl")
using .MMS


export SEM_Wave

mutable struct SEM_Wave
    #stores data for SEM approximation of the wave equation \frac{\partial^2 u}{\partial t^2} = \frac{\partial^2 u}{\partial x^2} - f(x)\cos(\omega t)
    #Dirichlet boundary conditions are used

    nodes::Vector{Float64}          #the endpoints of the elements
    PointsPerElement::Int64         #determines the degree of the Lagrange polynomials 
    x::Vector{Float64}              #the actual coordinates between nodes[1] and nodes[end] where we approximate our solution

    nsteps::Int64                   #number of steps we make in time
    Tend::Float64                   #the endpoint of our interval in time
    timestep::Float64               #size of the timesteps (= Tend/nsteps)
    
    uPrev::Vector{Float64}          #value of the displacement in the previous timepoint
    uNow::Vector{Float64}           #value of the displacement in the current timepoint
    uNext::Vector{Float64}          #value of the displacement in the next timepoint

    fVals::Vector{Float64}          #values of the driving term of the equation
    omega::Float64                  #frequency parameter omega
    
    QuadPoints::Vector{Float64}     #the points in [-1, 1] we use in our Gauss-Lobatto quadrature
    QuadWeights::Vector{Float64}    #the weights corresponding to the points in QuadPoints
    G::Array{Float64}               #a matrix required for the SEM method: G_{jm} = \sum_{l=0}^{N} w_l l'_m(\xi_l) l'_j(\xi_l); w_l \in QuadWeights, \xi_l \in QuadPoints, l_j's are the Lagrange polynomials
    inverseM::Vector{Float64}       #the inverse of the mass matrix, which is diagonal => stored as vector
    
    useMMS::Bool                    #determines whether an MMS test is run or not
    MMS_j::MMS_jet                  #stores the relevant data for the MMS test


    function SEM_Wave(
        nodes::Vector{Float64}, 
        nsteps::Int64,
        Tend::Float64,
        PointsPerElement::Int64,
        fVals::Vector{Float64}, 
        omega::Float64
        )
        
        QuadPoints, QuadWeights = gausslobatto(PointsPerElement)
        uPrev = zeros((length(nodes)-1)*(PointsPerElement-1) + 1)
        uNow = zeros((length(nodes)-1)*(PointsPerElement-1) + 1)
        uNext = zeros((length(nodes)-1)*(PointsPerElement-1) + 1)
        x = ConstructX(nodes, QuadPoints)

        useMMS = false

        timestep = Tend/nsteps

        #make G and the inverse of the mass matrix
        G = ConstructG(QuadPoints, QuadWeights)
        inverseM = ConstructInverseM(nodes, QuadWeights)

        #change some coefficients in the MMS
        MMS_j = MMS.MMS_jet(1, 2)
        MMS_j.coeff[1, 2] = 2.0
        MMS_j.coeff[2, 1] = exp(1)
        MMS_j.coeff[2, 3] = -2
        MMS_j.coeff[1, 3] = (1 + sqrt(5))/2

            new(nodes, PointsPerElement, x, nsteps, Tend, timestep,
                uPrev, uNow, uNext, fVals, omega, QuadPoints, QuadWeights, G, inverseM, useMMS, MMS_j)
        
    end

end


function MakeStep!(simul::SEM_Wave, stepnumber::Int64)
    

    simul.uNext = 2*simul.uNow - simul.uPrev + 
                    simul.timestep^2 * LaplaceCalculation(simul, simul.uNow) + 
                    simul.timestep^2 * ForcingTerm(simul, stepnumber)        
    
    ###########
    #enforce Dirichlet boundary conditions
    simul.uNext[1] = 0
    simul.uNext[end] = 0
    
    #=
    if simul.useMMS
        simul.uNext[1] = 
        simul.uNext[end] = 
    end 
    =#

    #update the variables
    simul.uPrev = simul.uNow
    simul.uNow = simul.uNext

end


function Simulate(nodes, nsteps, Tend, N, fVals, omega, uStart, uStartDer, animate = false, snapshotFrequency = 100, plotHeight = 10.0, animationName = "test.gif")
    #approximates a solution to the wave equation using SEM_Wave
    #creates a gif named animationName if animate == true
    #snapshotFrequency decides how many steps are made before a frame is saved in the gif
    #the y-axis shown in the gif is [-plotHeight, plotheight]

    #create a SEM_Wave with the specified parameter values
    simul = SEM_Wave(nodes, nsteps, Tend, N, fVals, omega)

    Initialise!(simul, uStart, uStartDer)
    

    #create an Animation if one is asked for
    if (animate == true)
       
        anim = Animation()

        plot(simul.x, simul.uNow, xlims=(nodes[1],nodes[end]), ylims=(-plotHeight, plotHeight))
        frame(anim)    
    
    end
    
    
    #make the steps
    for n=1:nsteps

        MakeStep!(simul, n)    

        #making the Animation is very slow compared to the actual calculations made here.
        #store a frame in our Animation once every snapshotFrequency steps
        if (animate == true && n%snapshotFrequency == 0)
            plot(simul.x, simul.uNow, xlims=(nodes[1],nodes[end]), ylims=(-plotHeight, plotHeight))
            frame(anim)
            #println(n)
        end

            
    end

    #finally make and store the gif
    if (animate == true)

        gif(anim, animationName, fps=40) #if you need the gif to run a different framerate, change the value of fps

    end

    return simul.uNow

end


function LaplaceCalculation(simul::SEM_Wave, u::Vector{Float64})

    #to = TimerOutput()
    N = simul.PointsPerElement     #number of quadrature points
    K = length(simul.nodes)-1      #number of elements
    degreesOfFreedom = K*(N-1) + 1
   
    laplaceVals = zeros(degreesOfFreedom)

    u_k = zeros(N)
    v_k = zeros(N)
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

        mul!(v_k, simul.G, u_k, 2/delta_x_k, 0)

        SetDegreesOfFreedom!(simul, k, laplaceVals, v_k, true)

    end

    #show(to)

    laplaceVals = -laplaceVals.*simul.inverseM      #the output matrix corresponds to -Laplace, so we need to change the sign

    return laplaceVals

end


function ForcingTerm(simul::SEM_Wave, stepnumber::Int64)

    N = simul.PointsPerElement     #number of quadrature points
    K = length(simul.nodes)-1      #number of elements
    degreesOfFreedom = K*(N-1) + 1
   
    forcingVals = zeros(degreesOfFreedom)

    F_k = zeros(N)
    vals_k = zeros(N)
    
    if simul.useMMS

        t = simul.timestep*stepnumber

        x_loc = zeros(N)

        z = (1 .+ simul.QuadPoints)

        for k = 1:K

            delta_x_k = simul.nodes[k+1] - simul.nodes[k]
            
            x_loc .= simul.nodes[k] .+ 0.5*z*delta_x_k

            for j = 1:N
                vals_k[j] = MMS.MMSfun(x_loc[j], t, 0, 2, simul.MMS_j) - MMS.MMSfun(x_loc[j], t, 2, 0, simul.MMS_j)
            end
            

            SetDegreesOfFreedom!(simul, k, forcingVals, vals_k, false)

        end

    else

        for k = 1:K

            delta_x_k = simul.nodes[k+1] - simul.nodes[k]

            F_k .= GetDegreesOfFreedom(simul, k, simul.fVals)

            vals_k = 0.5*delta_x_k * (simul.QuadWeights.*F_k)

            SetDegreesOfFreedom!(simul, k, forcingVals, vals_k, true)

        end

        forcingVals = forcingVals.*simul.inverseM * cos(simul.omega*(stepnumber-1)*simul.timestep)

    end

    return forcingVals

end


function Initialise!(simul::SEM_Wave, uStart::Vector{Float64}, uStartDer::Vector{Float64})
    
    simul.uNow = uStart

    simul.uPrev = uStart - simul.timestep*uStartDer + 0.5*simul.timestep^2*(LaplaceCalculation(simul, uStart) + ForcingTerm(simul, 1)) #stepnumber = 1 corresponds to t = 0

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


function ConstructInverseM(nodes::Vector{Float64}, QuadWeights::Vector{Float64})

    #constructs the inverse of the mass matrix

    K = length(nodes)-1
    N = length(QuadWeights)

    inverseM = zeros(K*(N-1) + 1)
    
    for k = 1:K
            
        start = (k-1)*(N-1)

        delta_x_k = nodes[k+1] - nodes[k]

        inverseM[(start + 1):(start + N)] += 0.5*QuadWeights*delta_x_k
        
    end

    inverseM = 1 ./inverseM

end


function GetDegreesOfFreedom(simul::SEM_Wave, k::Int64, u::Vector{Float64})
    
    start = (k-1)*(simul.PointsPerElement-1)
    return u[(start + 1):(start + simul.PointsPerElement)]

end


function SetDegreesOfFreedom!(simul::SEM_Wave, k::Int64, v::Vector{Float64}, v_k::Vector{Float64}, add::Bool)
    
    start = (k-1)*(simul.PointsPerElement-1)

    if add
        #v[(start + 1):(start + simul.PointsPerElement)] .+= v_k
        view(v, (start + 1):(start + simul.PointsPerElement)) .+= v_k
    else
        v[(start + 1):(start + simul.PointsPerElement)] .= v_k
    end

end

#this does not compare in simul.nodes[0] and simul.nodes[end].
#Change when the boundary conditions have been sorted out properly
function LaplaceMMS(simul::SEM_Wave, tol::Float64, type)   

    acceptableAccuracy = true

    if type == "polynomial"

        #make a polynomial which we will test against, and calculate its derivative symbolically
        #testPolynomial = Polynomial([1, 2, 3, 3, 5, 6, -7, 3, 1, -2, 1])
        #testPolynomial = Polynomial([0, 1, -1])^3 * Polynomial([0, 1])^4
        testPolynomial = Polynomial([0, 1, -1])^3 * Polynomial([0, 1])^4 *Polynomial([1, 2, -3, 4])
        #testPolynomial = Polynomial([0, 1])
        testPolynomial_xx = derivative(derivative(testPolynomial))

        testPolynomialValues = zeros(length(simul.uNow))
        testPolynomial_xxValues = zeros(length(testPolynomialValues))

        #evaluate the polynomial for the appropriate x-values
        for j = 1:length(simul.x)
            testPolynomialValues[j] = testPolynomial(simul.x[j])
            testPolynomial_xxValues[j] = testPolynomial_xx(simul.x[j])
        end

        #calculate the laplacian of this polynomial using LaplaceCalculation
        laplaceValues = LaplaceCalculation(simul, testPolynomialValues)

        #compare with the value from LaplaceCalculation, L^\infty-norm
        #might want to use the L^2-norm instead here
        diff = maximum(abs.(laplaceValues[2:end-1] - testPolynomial_xxValues[2 : end-1]))
        
        if (diff > tol)
            acceptableAccuracy = false
        end

        
    end

    #=
    if (acceptableAccuracy)
        println("LaplaceCalculation passed the MMS test")
    else
        println("LaplaceCalculation did not pass the MMS test")
    end
    =#

    return [laplaceValues[2:end-1], testPolynomial_xxValues[2:end-1]]
    
end


end