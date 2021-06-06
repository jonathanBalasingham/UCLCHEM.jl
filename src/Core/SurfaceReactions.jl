module SurfaceReactions
using .Constants

SILICATE_MU = 0.005D0 # Fraction of newly formed H2 that stays on the grain surface
SILICATE_E_S = 110.0D0 # Energy of the saddle point between a physisorbed and a chemisorbed site (K)
SILICATE_E_H2 = 320.0D0 # Desorption energy of H2 molecules (K)
SILICATE_E_HP = 450.0D0 # Desorption energy of physisorbed H atoms (K)
SILICATE_E_HC = 3.0D4   # Desorption energy of chemisorbed H atoms (K)
SILICATE_NU_H2 = 3.0D12 # Vibrational frequency of H2 molecules in surface sites (s^-1)
SILICATE_NU_HC = 1.3D13 # Vibrational frequency of H atoms in their surface sites (s^-1)

#Graphite grain properties for H2 Formation
GRAPHITE_MU = 0.005D0   # Fraction of newly formed H2 that stays on the grain surface
GRAPHITE_E_S = 260.0D0  # Energy of the saddle point between a physisorbed and a chemisorbed site (K)
GRAPHITE_E_H2 = 520.0D0 # Desorption energy of H2 molecules (K)
GRAPHITE_E_HP = 800.0D0 # Desorption energy of physisorbed H atoms (K)
GRAPHITE_E_HC = 3.0D4   # Desorption energy of chemisorbed H atoms (K)
GRAPHITE_NU_H2 = 3.0D12 # Vibrational frequency of H2 molecules in surface sites (s^-1)
GRAPHITE_NU_HC = 1.3D13 # Vibrational frequency of H atoms in their surface sites (s^-1)

GAS_DUST_MASS_RATIO = 100.0
GRAIN_RADIUS = 1.0d-5
GRAIN_DENSITY = 3.0 # Mass density of a dust grain
THERMAL_VEL = sqrt(8.0d0*K_BOLTZ/(PI*AMU)) #Thermal velocity without the factor of SQRT(T/m) where m is moelcular mass in amu

#reciprocal of fractional abundance of dust grains (we only divide by number density so better to store reciprocal)
GAS_DUST_DENSITY_RATIO = (4.0*PI*(GRAIN_RADIUS^3)*GRAIN_DENSITY*GAS_DUST_MASS_RATIO)/(3.0 * AMU)
#Grain area per h nuclei, calculated from average radius.
GRAIN_AREA_PER_H = 4.0*PI*GRAIN_RADIUS*GRAIN_RADIUS/GAS_DUST_DENSITY_RATIO

#Below are values for grain surface reactions
DIFFUSE_REACT_COMPETITION = true 
GRAINS_HAVE_ICE = true
CHEMICAL_BARRIER_THICKNESS = 1.40d-8  #gre Parameter used to compute the probability for a surface reaction with 
## activation energy to occur through quantum tunneling (Hasegawa et al. Eq 6 (1992).)
SURFACE_SITE_DENSITY = 1.5d15 # site density on one grain [cm-2]
VDIFF_PREFACTOR = 2.0*K_BOLTZ*SURFACE_SITE_DENSITY/PI/PI/AMU
NUM_SITES_PER_GRAIN = GRAIN_RADIUS*GRAIN_RADIUS*SURFACE_SITE_DENSITY*4.0*PI
end