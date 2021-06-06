using Flux
using UCLCHEM
using Sundials
using Dates

rfp = "./test/input/reactions.csv"
sfp = "./test/input/species.csv"
icfp = "./test/input/initcond0.csv"

T=10
zeta= 1.3e-17 
zeta=1.0
F_UV=1
A_v=2
omega = 0.5
E = 0.5
dens = 1e4
p = UCLCHEM.Parameters(zeta, omega, T, F_UV, A_v, E, dens)

#tspan = (0., 10^7 * 365. * 24. * 3600.)
tspan = (0., 10^4 * 365. * 24. * 3600.)

nw_prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan)
prob = ODEProblem(nw_prob.network, nw_prob.u0, tspan)
sol = solve(prob, CVODE_BDF(), saveat = 3600)

#data = transpose(reduce(hcat, sol.u))
x = sol.u[1:end-1]
y = sol.u[2:end]

chain = Flux.Chain(Dense(length(nw_prob.species),128), Dense(128, length(nw_prob.species)))

sum1_penalty = 1000

function predict_err(xk, yk)
    pred = chain(xk)
    Flux.Losses.mse(pred, yk) + sum1_penalty * min(0,abs(sum(pred) - 1))
end

loss() = sum(abs.(predict_err.(x,y)))

parameters = params(chain)
opt = ADAM(0.1)
data = Iterators.repeated((), 3000)

evalcb() = @show(loss())
@time Flux.train!(loss, parameters, data, opt, cb = Flux.throttle(evalcb, 1))

using BSON: @save
@save "roberson_network.bson" chain