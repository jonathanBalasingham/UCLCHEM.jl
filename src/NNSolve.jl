using DiffEqFlux
using NeuralPDE
using DifferentialEquations
using Flux

function nnsolve(prob::ChemicalNetworkProblem,
                 alg::NeuralPDE.NNODE;
                 tol::Float64=10^-20, 
                 maxiters::Int=20000, 
                 dt::Number = 20, 
                 dt_factor=1.4, 
                 time_factor=1.01, time_factor_pre_1000_years=10.)
 
    current_time = 0.0
    target_time = 1.1*time_factor_pre_1000_years
    u = Float64[]
    t = Float64[]
    year_1000 = 1000. * 3600. * 24. * 365.
    current_problem = ODEProblem{false}(prob.network, prob.u0, (current_time, target_time))

    while current_time <= prob.tspan[2]
        print(current_time)
        print(" ")
        sol = NeuralPDE.solve(current_problem, alg, dt=(target_time - current_time) / dt, 
                              verbose=true, abstol=tol, maxiters=maxiters, save_everystep=true)
        u = vcat(u, sol.u[2:end])
        t = vcat(t, sol.t[2:end])
        current_time = target_time
        if current_time > year_1000 target_time *= time_factor else target_time *= time_factor_pre_1000_years end
        current_problem = ODEProblem{false}(prob.network, sol.u[end], (current_time, target_time))
        dt *= dt_factor
    end
    ChemicalNetworkSolution(t,u, prob.species)
end


function nn_solve(prob::ChemicalNetworkProblem,
                alg::NeuralPDE.NNODE;
                tol::Float64=10^-20, 
                maxiters::Int=20000, 
                number_of_problems::Number = 20,
                dt_of_each_problem::Number = 1000,
                save_everystep=false)

    tspan = prob.tspan
    prob_tspan = (tspan[2] - tspan[1]) / number_of_problems
    tspans = tspan[1]:prob_tspan:tspan[2]
    u = Float64[]
    t = Float64[]
    current_problem = ODEProblem{false}(prob.network, prob.u0, (tspans[1], tspans[2]))

    for (i,j) in enumerate(tspans[2:end-1])
        println("Current time: ", j)
        sol = NeuralPDE.solve(current_problem, alg, dt=(tspans[i+1] - tspans[i]) / dt_of_each_problem, verbose=true, abstol=tol, maxiters=maxiters, save_everystep=true)
        current_problem = ODEProblem{false}(prob.network, sol.u[end], (tspans[i], tspans[i+1]))
        u = vcat(u, sol.u[2:end])
        t = vcat(t, sol.t[2:end])
    end
    ChemicalNetworkSolution(t,u, prob.species)
end