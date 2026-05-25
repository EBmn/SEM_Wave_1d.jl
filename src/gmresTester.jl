

using LinearAlgebra
using SparseArrays
using Plots
using FastGaussQuadrature
using Polynomials
using TimerOutputs
using IterativeSolvers

include("SEM_Wave_2d.jl")
using .SEM_Wave_2d


function main()

    N = 1000;

    A = I(N) + 0.05*rand(N, N)
    b = ones(N, 1)

    x = zeros(N, 1)
    #x, history = gmres!(x, A, b, verbose=true)

    x, history = gmres!(x, A, b, log=true, verbose = true);

    println(size(A))
    println(size(x))
    println(size(b))
    println(norm(A*x-b, 2))
    println(history[:resnorm])

end

