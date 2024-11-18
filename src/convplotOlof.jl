
using LinearAlgebra
using Plots
using FastGaussQuadrature
using LinearAlgebra
using TimerOutputs


include("SEM_Wave_1d.jl")

#=
    Produces convergence plots SEM_Wave for a fixed polynomial degree and increasing number of elements
    Max norm of differences used to illustrate the convergence.
=#

function main()

    # Problem parameters
    
    omega = 10*pi
    Tend = 0.9
    xl = 0.0
    xr = 1.0
    u0(x) = 0.0
    u0prim(x) = 0.0
    f(x) = exp(-300*(x-0.4)^2)-exp(-500*(x-0.75)^2)
    
    # Convergence parameters
    
    levels = 10            # Number of discretization levels
    NList = [3, 5, 6, 7]   # Polynomial orders
    M0 = 2;                # Smallest number of elements
    K0 = 30;               # Smallest number of timesteps    

    # Other initializations
    
    plt = plot()
    solprev = zeros(1)

    # Loop starts
    
    for N in NList
        
        diffs = zeros(levels-1)  # max norm difference between level j+1 and j
        Ms = zeros(levels-1)     # number of elements on level j+1 (matches with diffs)

        M = M0;   # Number of elements
	K = K0;   # Number of time steps
	
        for j = 1:levels

            # Grid for solver
	    
            nodes = collect(LinRange(xl, xr, M+1))
            quadpoints, ~ = gausslobatto(N)
            x = SEM_Wave_1d.ConstructX(nodes, quadpoints)
	    
            # Solve
	    
            sol = SEM_Wave_1d.Simulate(nodes, K, Tend, N, f.(x), omega, u0.(x), u0prim.(x), false, 0)

            if (j>1)
	       diffs[j-1] = maximum(abs.(sol[2*N-1:2*N-2:end-1]-solprev[N:N-1:end-1]))
	       Ms[j-1] = M
            end
 
            solprev = sol # Save previous solution to be able to compute diffs
	    M = M*2;
	    K = K*2;
        end

        # Plotting and post-processing

        plot!(Ms, diffs, xscale=:log10, yscale=:log10, xlims=(Ms[1], Ms[end]), label = "N = " * string(N), legend=:bottomleft)

        display(diffs[1:end-1]./diffs[2:end])  # Display ratios as a column
	
        println("done with N = " * string(N))        

    end

    #savefig(plt,"convplotNotAnalytic")
    display(plt)
end

