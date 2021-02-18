using Flux
using UCLCHEM
using Sundials
using Dates
using OrdinaryDiffEq
using ReservoirComputing
using ParameterizedFunctions

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
p = UCLCHEM.Parameters(zeta, omega, T, F_UV, A_v, E, dens)

#tspan = (0., 10^7 * 365. * 24. * 3600.)
tspan = (0., 10^0 * 365. * 24. * 3600.)

nw_prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan)
prob = ODEProblem(nw_prob.network, nw_prob.u0, tspan)
sol = solve(prob, CVODE_Adams())
#sol2 = UCLCHEM.solve(nw_prob, dt=10000)
println("Problem solved with CVODE")

v = sol.u
data = Matrix(hcat(v...))
shift = 1
train_len = 60000
predict_len = 10000
train = data[:, shift:shift+train_len-1]
test = data[:, shift+train_len:shift+train_len+predict_len-1]

approx_res_size = 300
radius = 2.1
degree = 5
activation = tanh
sigma = 0.7
beta = 0.0
alpha = 0.98 # leaking factor
nla_type = NLADefault()
extended_states = false


esn = ESN(approx_res_size,
    train,
    degree,
    radius,
    activation = activation, #default = tanh
    alpha = alpha, #default = 1.0
    sigma = sigma, #default = 0.1
    nla_type = nla_type, #default = NLADefault()
    extended_states = extended_states #default = false
    )

@time W_out = ESNtrain(esn, beta)
@time output = ESNpredict(esn, predict_len, W_out)


using Plots
plot(transpose(output[1:12,:]),layout=(4,3), label="predicted",size=(1200,800))
plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
