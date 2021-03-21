using Sundials, DiffEqBase, DifferentialEquations, Flux, DiffEqFlux, Plots
using UCLCHEM
using DiffEqCallbacks
using Dates

rfp = "./test/input/reactions_postNN.csv"
sfp = "./test/input/species_postNN.csv"
icfp = "./test/input/initcond1.csv"

T=10. 
zeta = 1.
omega = 0.5
F_UV=1.
A_v=2.
E = 0.5
dens = 1e4
params = UCLCHEM.Parameters(zeta, omega, T, F_UV, A_v, E, dens)

#tspan = (0., 10^7 * 365. * 24. * 3600.)
# tspan (0.1) takes ~ 30 to 50 seconds
# tspan (0,10) takes 60 to 80 seconds
# tspan (0,100) takes 5 to 8 minutes
tspan = (0., 1000.0)

p = UCLCHEM.formulate(sfp,rfp,icfp,params,tspan, rate_factor = 100)

u0 = p.u0

prob = ODEProblem(p.network,p.u0,tspan)
#sol = UCLCHEM.solve(p, solver=Rosenbrock23, maxiter=100000000000)
#sol = solve(prob,CVODE_Adams(),reltol=1e-12,abstol=1e-12)
#prob2 = ODEProblem(rober,sol[end],(100.0,0.0))
#sol = solve(prob,CVODE_Adams(),reltol=1e-12,abstol=1e-12)
#@show sol[end]-u0 #[-17.5445, -14.7706, 39.7985]
t = range(tspan[1],tspan[2],length=20)
ode_data = Array(solve(prob,CVODE_BDF(),saveat=t, reltol=1e-4,abstol=1e-20))

dudt = FastChain(FastDense(length(p.species),length(p.species),σ))
n_ode = NeuralODE(dudt,tspan,Euler(),dt=.01 ,saveat=t,reltol=1e-4,abstol=1e-20)
ps = Flux.params(n_ode)
pred = n_ode(prob.u0)
println("Plotting first scatter")
plot(t,transpose(ode_data),label="data")
plot!(t,transpose(pred),label="prediction")

function predict_n_ode()
    n_ode(u0)
end
loss_n_ode() = sum(abs2,ode_data .- predict_n_ode())

data = Iterators.repeated((), 500)
opt = ADAM(0.1)
println("Defining callback")

cb = function () #callback function to observe training
  display(loss_n_ode())
  # plot current prediction against data
  display(now())
  #pl = scatter(t,transpose(ode_data),label="data")
  #scatter!(pl,t,transpose(cur_pred),label="prediction")
  #display(plot(pl))
end
println("Running cb")
# Display the ODE with the initial parameter values.
cb()
println("training")
Flux.train!(loss_n_ode, ps, data, opt, cb=cb)

plot(t,transpose(ode_data),label="data")
plot(t,transpose(predict_n_ode()),label="pred")