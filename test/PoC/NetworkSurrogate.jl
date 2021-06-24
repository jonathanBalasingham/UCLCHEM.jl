using UCLCHEM, ReservoirComputing, Sundials, Surrogates

rfp = "./test/input/reactions_final.csv"
sfp = "./test/input/species_final.csv"
icfp = "./test/input/initcond_first_network.csv"

rfp = "./test/input/reactions.csv"
sfp = "./test/input/species.csv"
#icfp = "./test/input/initcond1.csv"


tspan = (0., 10^6 * 365. * 24. * 3600.)

rates_set_lower_bound = [0.5, 0.3, 5, 0.5, 0.5, 0.25, 1e4]
rates_set_upper_bound = [1.5, 1.0, 15, 1.5, 1.5, .75, 1e4]
T=10 
CR= 1.3e-17 
FUV=1
Av=10

T=10
zeta= 1.3e-17 
#zeta=1.0
F_UV=1
A_v=1.0002E+01
omega = 0.6
E = 0.5
dens = 1e4
pa = Parameters(zeta, omega, T, F_UV, A_v, E, dens)

params_from_first = [1.,1.,10.,1.,10, 1.3e-17,1e4]
struct Parameters
    zeta::Float64
    omega::Float64
    T::Float64
    F_UV::Float64
    A_v::Float64
    E::Float64
    density::Float64
end
x = sample(30, rates_set_lower_bound, rates_set_upper_bound, SobolSample())



#= ESN Settings =#
res_size = 1000
radius = .6
degree = 600
activation = tanh
alpha = .6
sigma = .1
nla_type = NLADefault()
extended_states = false
beta = 0.000001

y = []

i = 1
parameters = x[begin]
pa = Parameters(parameters...)
#p = UCLCHEM.formulate(sfp,rfp,icfp,pa,tspan, rate_factor = 1)
p = formulate_all(rfp, icfp, pa)
sol = solve(p, solver=CVODE_BDF)

train_subset = hcat(sol.u...) |> Matrix
#train_subset[train_subset .<= 0.0] .= 1e-60

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
successful_indx = zeros(length(x))

for parameters in x
    try
        println(i)
        pa = Parameters(parameters...)
        p = formulate_all(rfp, icfp, pa)
        #p = UCLCHEM.formulate(sfp,rfp,icfp,pa,tspan, rate_factor = 1)
        @time sol = solve(p, solver=CVODE_BDF);   
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
        successful_indx[i] = i   
    catch
        println("Failed with parameters: $pa")
    end  
end

radial_surrogate = RadialBasis(x, y, rates_set_lower_bound, rates_set_upper_bound)

test_parameters = [0.6, 0.7, 9, 0.6, 0.9, 0.6, 1.e4]
test_W_out = reshape(radial_surrogate(test_parameters), length(p.species), :)

pa = Parameters(parameters...)
#p = UCLCHEM.formulate(sfp,rfp,icfp,pa,tspan, rate_factor = 1)
p = formulate_all(rfp, icfp, pa)
sol = solve(p, solver=CVODE_BDF)
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

#plot(output[1:20, :]',layout=(4,5), label="Predicted", legend=:outertopright, size=(1800, 800))
#title!(string(test_parameters))

#plot!(test_output[1:20, :]',layout=(4,5), label="Interpolated", size=(1800, 1200))

#plot!(train_subset[1:20, 2:end]', layout=(4,5), label="Ground Truth", size=(1800, 1200))
#xaxis!(:log10)

#output[output .<= 0.0] .= 1e-60
#train_subset[train_subset .<= 0.0] .= 1e-60
#test_output[test_output .<= 0.0] .= 1e-60

l = @layout [a{0.01h}; grid(2,2)]
#[plots[i+1] = plot(sol.t,train_subset[i, :],framestyle=:default,label="ground truth",title=p.species[i]) for i in 1:length(p.species)]
plots = [plot(sol.t ./ (3600 * 24 * 365),train_subset[i, :],label="ground truth",title=p.species[i], size=(300,300), yaxis=:log10, xaxis=:log10) for i in 1:length(p.species)]
plot(plots..., size=(1200,1000))


using Plots
for (i, species) in enumerate(p.species)
    plot(sol.t ./ (3600 * 24 * 365), output[i, :], title=species, label="Predicted", legend=:outertopright)
    plot!(sol.t ./ (3600 * 24 * 365), train_subset[i, :], title=species, label="Groud Truth", legend=:outertopright)
    plot!(sol.t ./ (3600 * 24 * 365), test_output[i, :], title=species, label="Interpolated", legend=:outertopright)
    xticks!([10^0, 10,10^2,10^3,10^4,10^5,10^6])
    xaxis!(:log10)
    ylims!((10^-30, 1))
    yaxis!(:log10)
    savefig("./output/$species.png")
end
