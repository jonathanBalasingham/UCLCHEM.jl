using DiffEqFlux
using NeuralPDE
using DifferentialEquations
using Flux
include("NNODE.jl")

function nnsolve(prob::ODEProblem, chain=Flux.Chain(Dense(1, 5, σ), Dense(5, 1)), tol::Float64=1e-10, maxiters::Int=1000, dt::Number = 1000)
    opt = Flux.ADAM(0.1, (0.9, 0.95))
    sol2 = NeuralPDE.solve(prob, NNODE(chain, opt), dt=dt, 
                       verbose=true,
                       abstol=tol, maxiters=maxiters)
    sol2
end