using UCLCHEM

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
UCLCHEM.visualize(sol, species=["H","H2", "CO", "CH3OH", "#CO", "#CH3OH"])

sol2 = UCLCHEM.nnsolve(prob)

