using Sundials
using DifferentialEquations

function solve(prob::ChemicalNetworkProblem, tol::Float64=1e-10, maxiter::Int=10000, dt::Number = 1000)
    current_time = 0.0
    target_time = 1.1
    u = Float64[]
    t = Float64[]
    current_problem = ODEProblem{false}(prob.network, prob.u0, (current_time, target_time))
    while current_time <= prob.tspan[2]
        print(current_time)
        print(" ")
        sol = DifferentialEquations.solve(current_problem, CVODE_Adams(), maxiter=maxiter)
        u = vcat(u, sol.u[2:end])
        t = vcat(t, sol.t[2:end])
        current_time = target_time
        target_time *= 1.1
        current_problem = ODEProblem{false}(prob.network, sol.u[end], (current_time, target_time))
    end
    ChemicalNetworkSolution(t,u, prob.species)
end