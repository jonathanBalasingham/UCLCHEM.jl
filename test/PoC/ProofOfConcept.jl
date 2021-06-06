using DifferentialEquations
using Sundials
using Plots
using ReservoirComputing

function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁+k₃*y₂*y₃
  du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  du[3] =  k₂*y₂^2
  nothing
end

tspan = (0.0,1e5)
prob = ODEProblem(rober,[0.5,0.3,0.2],tspan,[0.01,3e3,1e7])
sol = solve(prob, CVODE_BDF()) #, saveat=sol.t)
#plot(sol,tspan=(1e-2,1e5),xscale=:log10)

sol_as_matrix = hcat(sol.u...) |> Matrix


res_size = 500
radius = .2
degree = 100
activation = tanh
alpha = 1.
sigma = .1
nla_type = NLADefault()
extended_states = false
beta = 0.000001

#= Train/Test interweaved  =#
train_subset = sol_as_matrix[:, 1:2:end-1]
test_subset = sol_as_matrix[:, 2:2:end]

#= Train/Test sequential  =#
#train_subset = sol_as_matrix[:, 1:100]
#test_subset = sol_as_matrix[:, 101:end]

#= Log Abundances =#
#train_subset = log10.(train_subset)
#test_subset = log10.(test_subset)

esn = ESN(res_size,
          train_subset,
          degree,
          radius,
          activation = activation, #default = tanh
          alpha = alpha, #default = 1.0
          sigma = sigma, #default = 0.1
          nla_type = nla_type, #default = NLADefault()
          extended_states = extended_states)


@time W_out = ESNtrain(esn, beta)
# reset the states
#esn.states[:, end] = esn.states[:, 2]
fake_state = zeros(Float64, esn.res_size, 1)
y = test_subset[:, 1]
esn.states[:,end] =  (1-alpha).*fake_state + alpha*activation.((esn.W*fake_state)+(esn.W_in*y))
@time output = ESNpredict(esn, size(test_subset, 2)-1, W_out)

plot(output',layout=(3,1), label="predicted")
#plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
plot!(test_subset[:, 2:end]', layout=(3,1), label="actual")
xaxis!(:log10)

#= Testing the Weight prediction =# 

chain = 0
BSON.@load "./test/PoC/roberson_network.bson" chain

rates = [0.65, 3e2, 763112]
prob = ODEProblem(rober, [0.5,0.3,0.2],tspan, rates)
sol = solve(prob, CVODE_BDF()) #, saveat=sol.t)

pred_W_out = chain(rates)
pred_W_out = reshape(pred_W_out, (498, 3))' |> Matrix .|> Float64
sol_as_matrix = hcat(sol.u...) |> Matrix 


train_subset = sol_as_matrix[:, 1:2:end-1]
test_subset = sol_as_matrix[:, 2:2:end]

esn = ESN(res_size,
          train_subset,
          degree,
          radius,
          activation = activation, #default = tanh
          alpha = alpha, #default = 1.0
          sigma = sigma, #default = 0.1
          nla_type = nla_type, #default = NLADefault()
          extended_states = extended_states)


fake_state = zeros(Float64, esn.res_size, 1)
y = test_subset[:, 1]
esn.states[:,end] =  (1-alpha).*fake_state + alpha*activation.((esn.W*fake_state)+(esn.W_in*y))

@time output = ESNpredict(esn, size(test_subset, 2)-1, pred_W_out)

plot(output',layout=(3,1), label="predicted")
#plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
plot!(test_subset[:, 2:end]', layout=(3,1), label="actual")
xaxis!(:log10)