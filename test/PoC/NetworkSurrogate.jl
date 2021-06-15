using UCLCHEM, ReservoirComputing, Sundials, Surrogates

rfp = "./test/input/reactions_final.csv"
sfp = "./test/input/species_final.csv"
icfp = "./test/input/initcond0.csv"

#rfp = "./test/input/reactions.csv"
#sfp = "./test/input/species.csv"
#icfp = "./test/input/initcond1.csv"


tspan = (0., 10^6 * 365. * 24. * 3600.)

rates_set_lower_bound = [0.5, 0.3, 5, 0.5, 0.5, 0.25, 1e4]
rates_set_upper_bound = [1.5, 1.0, 15, 1.5, 1.5, .75, 1e5]

x = sample(500, rates_set_lower_bound, rates_set_upper_bound, SobolSample())



#= ESN Settings =#
res_size = 1000
radius = .6
degree = 600
activation = tanh
alpha = 1.
sigma = .1
nla_type = NLADefault()
extended_states = false
beta = 0.000001

y = []

i = 1
parameters = x[begin]
pa = UCLCHEM.Parameters(parameters...)
p = UCLCHEM.formulate(sfp,rfp,icfp,pa,tspan, rate_factor = 1)
sol = UCLCHEM.solve(p, solver=QNDF1)

train_subset = hcat(sol.u...) |> Matrix
esn = ESN(res_size,
        train_subset,
        degree,
        radius,
        activation = activation,
        alpha = alpha, 
        sigma = sigma, 
        nla_type = nla_type, 
        extended_states = extended_states)

W = esn.W
W_in = esn.W_in

for parameters in x[1:50]   
    println(i)
    pa = UCLCHEM.Parameters(parameters...)
    p = UCLCHEM.formulate(sfp,rfp,icfp,pa,tspan, rate_factor = 1)
    @time sol = UCLCHEM.solve(p, solver=QNDF1);   
    train_subset = hcat(sol.u...) |> Matrix

    esn = ESN(W,
            train_subset,
            W_in,
            activation = activation,
            alpha = alpha, 
            nla_type = nla_type, 
            extended_states = extended_states)

    @time W_out = ESNtrain(esn, beta)
    flattened_W_out = reshape(W_out, :, 1)
    push!(y, flattened_W_out)
    i += 1
end

radial_surrogate = RadialBasis(x[1:50], y, rates_set_lower_bound, rates_set_upper_bound)

test_parameters = [0.6, 0.2, 9, 0.6, 0.9, 0.45, 1.2e4]
test_W_out = reshape(radial_surrogate(test_parameters), 23, :)

pa = UCLCHEM.Parameters(parameters...)
p = UCLCHEM.formulate(sfp,rfp,icfp,pa,tspan, rate_factor = 1)
sol = UCLCHEM.solve(p, solver=QNDF1)
train_subset = hcat(sol.u...) |> Matrix

esn = ESN(W,
train_subset,
W_in,
activation = activation,
alpha = alpha, 
nla_type = nla_type, 
extended_states = extended_states)

@time W_out = ESNtrain(esn, beta)

@time test_output = ESNfitted(esn, test_W_out)
@time output = ESNfitted(esn, W_out)

plot(output[1:20, :]',layout=(4,5), label="Predicted", legend=:outertopright, size=(1800, 800))
#title!(string(test_parameters))

plot!(test_output[1:20, :]',layout=(4,5), label="Interpolated", size=(1800, 1200))

plot!(train_subset[1:20, 2:end]', layout=(4,5), label="Ground Truth", size=(1800, 1200))
xaxis!(:log10)

using Plots
for (i, species) in enumerate(p.species)
    plot(sol.t ./ (3600 * 24 * 365), output[i, :], title=species, label="Predicted", legend=:outertopright)
    plot!(sol.t ./ (3600 * 24 * 365), train_subset[i, :], title=species, label="Groud Truth", legend=:outertopright)
    xticks!([10^0, 10,10^2,10^3,10^4,10^5,10^6])
    xaxis!(:log10)
    savefig("./output/$species.png")
end
