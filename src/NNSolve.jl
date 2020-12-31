using DiffEqFlux
using NeuralPDE
using DifferentialEquations
using Flux
include("NNODE.jl")

function nnsolve(prob::ChemicalNetworkProblem; chain=Flux.Chain(Dense(1, 5, σ), Dense(5, 1)), tol::Float64=1e-10, maxiters::Int=10000, dt::Number = 1000)
 
    opt = Flux.ADAM(0.1, (0.9, 0.95))
    current_time = 0.0
    target_time = 1.1
    u = Float64[]
    t = Float64[]
    year_1000 = 1000. * 3600. * 24. * 365.
    current_problem = ODEProblem{true}(prob.network, prob.u0, (current_time, target_time))

    while current_time <= prob.tspan[2]
        print(current_time)
        print(" ")
        sol = NeuralPDE.solve(current_problem, NNODE(chain, opt), dt=dt, 
        verbose=true,
        abstol=tol, maxiters=maxiters)
        u = vcat(u, sol.u[2:end])
        t = vcat(t, sol.t[2:end])
        current_time = target_time
        if current_time > year_1000 target_time *= 1.1 else target_time *= 10. end
        current_problem = ODEProblem{true}(prob.network, sol.u[end], (current_time, target_time))
    end
    ChemicalNetworkSolution(t,u, prob.species)
end