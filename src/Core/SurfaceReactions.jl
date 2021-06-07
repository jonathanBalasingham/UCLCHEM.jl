module SurfaceReactions

include("Constants.jl")
using .Constants

SILICATE_MU = 0.005e0 # Fraction of newly formed H2 that stays on the grain surface
SILICATE_E_S = 110.0e0 # Energy of the saddle point between a physisorbed and a chemisorbed site (K)
SILICATE_E_H2 = 320.0e0 # Desorption energy of H2 molecules (K)
SILICATE_E_HP = 450.0e0 # Desorption energy of physisorbed H atoms (K)
SILICATE_E_HC = 3.0e4   # Desorption energy of chemisorbed H atoms (K)
SILICATE_NU_H2 = 3.0e12 # Vibrational frequency of H2 molecules in surface sites (s^-1)
SILICATE_NU_HC = 1.3e13 # Vibrational frequency of H atoms in their surface sites (s^-1)

#Graphite grain properties for H2 Formation
GRAPHITE_MU = 0.005e0   # Fraction of newly formed H2 that stays on the grain surface
GRAPHITE_E_S = 260.0e0  # Energy of the saddle point between a physisorbed and a chemisorbed site (K)
GRAPHITE_E_H2 = 520.0e0 # Desorption energy of H2 molecules (K)
GRAPHITE_E_HP = 800.0e0 # Desorption energy of physisorbed H atoms (K)
GRAPHITE_E_HC = 3.0e4   # Desorption energy of chemisorbed H atoms (K)
GRAPHITE_NU_H2 = 3.0e12 # Vibrational frequency of H2 molecules in surface sites (s^-1)
GRAPHITE_NU_HC = 1.3e13 # Vibrational frequency of H atoms in their surface sites (s^-1)

GAS_DUST_MASS_RATIO = 100.0
GRAIN_RADIUS = 1.0e-5
GRAIN_DENSITY = 3.0 # Mass density of a dust grain
THERMAL_VEL = sqrt(8.0 * K_BOLTZ/(pi*AMU)) #Thermal velocity without the factor of SQRT(T/m) where m is moelcular mass in amu

#reciprocal of fractional abundance of dust grains (we only divide by number density so better to store reciprocal)
GAS_DUST_DENSITY_RATIO = (4.0*pi*(GRAIN_RADIUS^3)*GRAIN_DENSITY*GAS_DUST_MASS_RATIO)/(3.0 * AMU)
#Grain area per h nuclei, calculated from average radius.
GRAIN_AREA_PER_H = 4.0*pi*GRAIN_RADIUS*GRAIN_RADIUS/GAS_DUST_DENSITY_RATIO

#Below are values for grain surface reactions
DIFFUSE_REACT_COMPETITION = true 
GRAINS_HAVE_ICE = true
CHEMICAL_BARRIER_THICKNESS = 1.40e-8  #gre Parameter used to compute the probability for a surface reaction with 
## activation energy to occur through quantum tunneling (Hasegawa et al. Eq 6 (1992).)
SURFACE_SITE_DENSITY = 1.5e15 # site density on one grain [cm-2]
VDIFF_PREFACTOR = 2.0*K_BOLTZ*SURFACE_SITE_DENSITY/pi/pi/AMU
NUM_SITES_PER_GRAIN = GRAIN_RADIUS*GRAIN_RADIUS*SURFACE_SITE_DENSITY*4.0*pi
"""
h2FormRate
gasTemp
dustTemp
"""
function h2FormRate(gasTemp, dustTemp)
    #  Mean thermal velocity of hydrogen atoms (cm s^-1)
    thermalVel = 1.45e5 * sqrt(gasTemp / 1.0e2)

    #  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
    #  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
    stickingCoefficient = 1.0 / (1.0 + 0.04 * sqrt(gasTemp + dustTemp) + 0.2 * (gasTemp / 1.0e2)+ 0.08 * (gasTemp / 1.0e2)^2)
    FLUX = 1.0e-10 # Flux of H atoms in monolayers per second (mLy s^-1)
    #Our cross-sectional area per H is different to Cazaux and Tielens so we scale theirs
    #by the fraction of our standard value. 
    CROSS_SECTION_SCALE = GRAIN_AREA_PER_H / 1.660539e-21

    SILICATE_CROSS_SECTION = 8.473e-22 * CROSS_SECTION_SCALE # Silicate grain cross section per H nucleus (cm^-2/nucleus)
    GRAPHITE_CROSS_SECTION = 7.908e-22 * CROSS_SECTION_SCALE # Graphite grain cross section per H nucleus (cm^-2/nucleus)

    FACTOR1 = SILICATE_MU * FLUX / (2*SILICATE_NU_H2 * exp(-SILICATE_E_H2/dustTemp))
    FACTOR2 = 1.0*(1.0+sqrt((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))^2 /4.0*exp(-SILICATE_E_S/dustTemp)
    EPSILON = 1.0/(1.0+SILICATE_NU_HC/(2*FLUX) * exp(-1.5*SILICATE_E_HC/dustTemp) * (1.0+sqrt((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))^2)
    silicateFormationEfficiency = 1.0 / (1.0 + FACTOR1 + FACTOR2) * EPSILON

    FACTOR1 = GRAPHITE_MU*FLUX/(2*GRAPHITE_NU_H2 * exp(-GRAPHITE_E_H2/dustTemp))
    FACTOR2 = 1.0*(1.0+sqrt((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))^2 / 4.0*exp(-GRAPHITE_E_S/dustTemp)
    EPSILON = 1.0/(1.0+GRAPHITE_NU_HC/(2*FLUX) * exp(-1.5*GRAPHITE_E_HC/dustTemp) * (1.0+sqrt((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))^2)
    graphiteFormationEfficiency = 1.0/(1.0+FACTOR1+FACTOR2)*EPSILON

    #Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
    #Cazaux & Tielens (2004, ApJ, 604, 222)
    0.5*thermalVel*(SILICATE_CROSS_SECTION*silicateFormationEfficiency + GRAPHITE_CROSS_SECTION*graphiteFormationEfficiency)*stickingCoefficient
end

end