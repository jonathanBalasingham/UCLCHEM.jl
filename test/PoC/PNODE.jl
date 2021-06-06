using DiffEqFlux
using Sundials
using Flux
using DifferentialEquations


function rober(du,u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    du[1] = -k₁*y₁+k₃*y₂*y₃
    du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
    du[3] =  k₂*y₂^2
    nothing
end

number_of_species = 3
number_of_reactions = 3
total_inputs = number_of_reactions + number_of_species

parameters = [0.01*1.0,3e3*1.0,1e7*1.0]
u0 = [0.5,0.3,0.2]
tspan = (0.0, 1e2)
prob = ODEProblem(rober, u0,tspan, parameters)
sol = solve(prob, CVODE_BDF())

ode_data = Array(sol)
rates = repeat(parameters', size(ode_data, 2))' |> Matrix
surrogate_data = vcat(ode_data, rates)

chain = Chain(Dense(total_inputs,50),
             Dense(50,total_inputs))

n_ode = NeuralODE(chain, tspan, Tsit5(), saveat=sol.t, reltol=1e-4,abstol=1e-20)
ps = Flux.params(n_ode)

function predict_n_ode()
    n_ode([u0; parameters])
end
n = length(surrogate_data)
loss_n_ode() = sum(abs2, surrogate_data .- predict_n_ode()) / n

data = Iterators.repeated((), 500)
opt = ADAM(0.5)

callback() = show(loss_n_ode())

using Flux: @epochs
@epochs 50 Flux.train!(loss_n_ode, ps, data, opt, cb=callback)

