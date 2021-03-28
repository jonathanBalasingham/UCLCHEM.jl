using Sundials
using DifferentialEquations
using DelimitedFiles

function solve(prob::ChemicalNetworkProblem; 
                abstol::Float64=10^-20, 
                reltol=10^-4, 
                maxiter::Int=10000, 
                dt::Number = 100, 
                solver=CVODE_BDF,
                time_factor=1.1, time_factor_pre_1000_years=10.)

    current_time = 0.0
    target_time = 1.1
    u = Float64[]
    t = Float64[]
    year_1000 = 1000. * 3600. * 24. * 365.
    current_problem = ODEProblem{true}(prob.network, prob.u0, (current_time, target_time))
    while current_time <= prob.tspan[2]
        print(current_time)
        print(" ")
        sol = DifferentialEquations.solve(current_problem, solver(), maxiter=maxiter, reltol=reltol, abstol=abstol)
        u = vcat(u, sol.u[2:end])
        t = vcat(t, sol.t[2:end])
        current_time = target_time
        if current_time > year_1000
             target_time *= time_factor 
        else 
            target_time *= time_factor_pre_1000_years 
        end
        current_problem = remake(current_problem, tspan=(current_time, target_time), u0=sol.u[end])
    end
    ChemicalNetworkSolution(t,u, prob.species)
end

function solve(prob::ChemicalNetworkProblem,
               filepath::AbstractString; 
               abstol::Float64=10^-20, 
               reltol=10^-4, 
               maxiter::Int=10000,
               saveat=60, 
               solver=CVODE_BDF)

    current_time = 0.0
    target_time = 1.1
    year_1000 = 1000. * 3600. * 24. * 365.
    time_factor_pre_1000_years = 10.
    time_factor_post_1000_years = 1.1
    current_problem = ODEProblem(prob.network, prob.u0, (current_time, target_time))

    while current_time <= prob.tspan[2]
        println(current_time)
        sol = DifferentialEquations.solve(current_problem, solver(), maxiter=maxiter, reltol=reltol, abstol=abstol)

        data = Matrix(hcat(sol.u...))

        open(filepath, "a") do io
            writedlm(io, [sol.t data'], ',')
        end

        current_time = target_time
        if current_time > year_1000
            target_time *= time_factor_post_1000_years 
        else 
            target_time *= time_factor_pre_1000_years 
        end
        current_problem = remake(current_problem, tspan=(current_time, target_time), u0=sol.u[end])
    end

end