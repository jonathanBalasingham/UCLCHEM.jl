using UCLCHEM
using Plots
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

tspan = (0., 10^7 * 365. * 24. * 3600.)
prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan)
sol = UCLCHEM.solve(prob)
UCLCHEM.visualize(sol, species=["H","H2", "H2O", "NH3", "CO", "#CO", "#CH3OH"])

savefig(UCLCHEM.visualize(sol, species=["H","H2", "CO", "#CO", "#CH3OH"]), "CVODE_solution_10k.png")

sol2 = UCLCHEM.nnsolve(prob, maxiters=10000)
UCLCHEM.visualize(sol2, species=["H","H2", "CO", "#CO", "#CH3OH"])


sol2 = UCLCHEM.nnsolve(prob, maxiters=1000)
savefig(UCLCHEM.visualize(sol2, species=["H","H2", "CO", "#CO", "#CH3OH"]), "NN_solution_1K.png")
sol2 = UCLCHEM.nnsolve(prob, maxiters=10000)
savefig(UCLCHEM.visualize(sol2, species=["H","H2", "CO", "#CO", "#CH3OH"]), "NN_solution_10K.png")
sol2 = UCLCHEM.nnsolve(prob, maxiters=100000)
savefig(UCLCHEM.visualize(sol2, species=["H","H2", "CO", "#CO", "#CH3OH"]), "NN_solution_100K.png")

sol2 = UCLCHEM.nnsolve(prob, maxiters=100000000, dt_factor=1.8)
savefig(UCLCHEM.visualize(sol2, species=["H","H2", "CO", "#CO", "#CH3OH"]), "NN_solution_100M.png")


## Old Network
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
tspan = (0., 10^7 * 365. * 24. * 3600. / 100000)

prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan, rate_factor = 100000)

using Flux, DiffEqFlux
interior_nodes = 10
chain = DiffEqFlux.FastChain(DiffEqFlux.FastDense(1, interior_nodes, σ), DiffEqFlux.FastDense(interior_nodes, length(prob.species)))
opt = Flux.ADAM(0.01, (0.9, 0.95))
using NeuralPDE
sol2 = UCLCHEM.nnsolve(prob, NNODE(chain, opt))

sol = UCLCHEM.solve(prob, dt=100)


using Plots
UCLCHEM.visualize(sol, species=["H","H2", "CO", "CO2", "E-", "H2O"])
savefig(UCLCHEM.visualize(sol, species=["H","H2", "CO", "CO2", "E-", "H2O"]), "CVODE_solution_10k_first_network.png")

# these are tests to see if we can obtain convergence in 10,000 iterations by 
# decreasing the time factor

# umist
rfp = "./test/input/umist12.csv"
sfp = "./test/input/default_species.csv"
icfp = "./test/input/initcond1.csv"

T=10. 
zeta = 1.
omega = 0.5
F_UV=1.
A_v=2.
E = 0.5
dens = 1e4
p = UCLCHEM.Parameters(zeta, omega, T, F_UV, A_v, E, dens)

tspan = (0., 10^7 * 365. * 24. * 3600.)
prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan)
sol = UCLCHEM.solve(prob)

