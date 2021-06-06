using Sundials
using DifferentialEquations
<<<<<<< HEAD
using DelimitedFiles
using ProgressMeter
||||||| merged common ancestors
=======
using DelimitedFiles
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77

function solve(prob::ChemicalNetworkProblem; 
                abstol::Float64=10^-20, 
                reltol=10^-4, 
                maxiter::Int=10000, 
<<<<<<< HEAD
                solver=CVODE_BDF,
||||||| merged common ancestors
                dt::Number = 100, 
                solver=CVODE_Adams,
=======
                dt::Number = 100, 
                solver=CVODE_BDF,
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77
                time_factor=1.1, time_factor_pre_1000_years=10.)

    current_time = 0.0
    target_time = 1.1
    u = Float64[]
    t = Float64[]
    year_1000 = 1000. * 3600. * 24. * 365.
    current_problem = ODEProblem{true}(prob.network, prob.u0, (current_time, target_time))

    @async p = ProgressThresh(prob.tspan[2], 0)

    while current_time <= prob.tspan[2]
<<<<<<< HEAD
        sol = DifferentialEquations.solve(current_problem, solver(), maxiter=maxiter, reltol=reltol, abstol=abstol)
||||||| merged common ancestors
        print(current_time)
        print(" ")
        sol = DifferentialEquations.solve(current_problem, solver(), maxiter=maxiter, saveat=(target_time - current_time)/dt, reltol=reltol, abstol=abstol)
=======
        print(current_time)
        print(" ")
        sol = DifferentialEquations.solve(current_problem, solver(), maxiter=maxiter, reltol=reltol, abstol=abstol)
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77
        u = vcat(u, sol.u[2:end])
        t = vcat(t, sol.t[2:end])
        current_time = target_time
<<<<<<< HEAD
        if current_time > year_1000
             target_time *= time_factor 
        else 
            target_time *= time_factor_pre_1000_years 
        end
        @async update!(p, current_time |> floor |> Integer)
||||||| merged common ancestors
        if current_time > year_1000 target_time *= time_factor else target_time *= 10. end
        #current_problem = ODEProblem{true}(prob.network, sol.u[end], (current_time, target_time))
=======
        if current_time > year_1000
             target_time *= time_factor 
        else 
            target_time *= time_factor_pre_1000_years 
        end
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77
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
    year_in_secs = 3600 * 24 * 365

    open(filepath, "a") do io
        writedlm(io, Array([0.0; prob.u0]'), ',')
    end

    while current_time <= prob.tspan[2]
        println(current_time)
        

        if current_time > year_in_secs
            sol = DifferentialEquations.solve(current_problem, solver(), maxiter=maxiter, reltol=reltol, abstol=abstol, saveat=year_in_secs, dt=1.)
            data = Matrix(hcat(sol.u...))

            open(filepath, "a") do io
                writedlm(io, [sol.t data'], ',')
            end
        else
            sol = DifferentialEquations.solve(current_problem, solver(), maxiter=maxiter, reltol=reltol, abstol=abstol, dt=1.)
            println(length(sol))
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