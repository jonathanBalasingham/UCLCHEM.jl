module Rates

ubound(array :: Array, dim) = last(axes(array)[dim])
lbound(array :: Array, dim) = first(axes(array)[dim]

function calculateReactionRates()
    idx1=crpReacs[1]
    idx2=crpReacs[2]
    rate[idx1:idx2]=alpha[idx1:idx2]*zeta

    # UV photons, radfield has (factor of 1.7 conversion from habing to Draine)
    idx1=photonReacs[1]
    idx2=photonReacs[2]
    rate[idx1:idx2] = alpha[idx1:idx2]*exp(-gama[idx1:idx2]*av[dstep])*radfield/1.7

    # Reactions involving cosmic ray induced photon
    idx1=crphotReacs[1]
    idx2=crphotReacs[2]
    rate[idx1:idx2] .= alpha[idx1:idx2].*gama[idx1:idx2].*1.0/(1.0-omega).*zeta.*(gasTemp[dstep]/300).^beta[idx1:idx2]


    # freeze out only happens if fr>0 and depending on evap choice 
    idx1=freezeReacs[1]
    idx2=freezeReacs[2]
    cion=1.0+16.71d-4/(GRAIN_RADIUS*gasTemp[dstep])
    if fr == 0.0 || gasTemp[dstep] > 30.0
        rate[idx1:idx2]=0.0
    else
        for j in idx1:idx2

            if re1[j] == nelec
                rate[j] = THERMAL_VEL*alpha[j]*GRAIN_AREA_PER_H*fr*cion
            else
                # taken from Rawlings et al. 1992
                # 0.25 factor converts surface area of a spherical grain to cross-sectional area
                # Saves recalculating pi.r^2/n_H when we have 4pi.r^2/n_H for other parts of the code
                rate[j] = alpha[j]*THERMAL_VEL*sqrt(gasTemp[dstep]/mass[re1[j]])*0.25*GRAIN_AREA_PER_H*fr
                if beta[j] # = 0.0 rate[j] = rate[j]*cion end
            end
        end
    end

    # Desorption due to energy released by H2 Formations
    idx1=desoh2Reacs[1]
    idx2=desoh2Reacs[2]
    if desorb == 1 && h2desorb == 1 && gama[j] <= ebmaxh2 && mantle[dstep] >= 1.0e-30
        # Epsilon is efficieny of this process, number of molecules removed per event
        # h2form is formation rate of h2, dependent on hydrogen abundance. 
        rate[idx1:idx2] = epsilon*h2form*abund[nh,dstep]*1.0 / mantle[dstep]
    else
        rate[idx1:idx2] = 0.0
    end
   
    # Desorption due to energy from cosmic rays
    idx1=descrReacs[1]
    idx2=descrReacs[2]
    if desorb == 1 && crdesorb == 1 && mantle[dstep] >= 1e-30 && gama[j] <= ebmaxcr 
        # 4*pi*zeta = total CR flux. 1.64d-4 is iron to proton ratio of CR
        # as iron nuclei are main cause of CR heating.
        # GRAIN_AREA_PER_H is the total area per hydrogen atom. ie total grain area per cubic cm when multiplied by density.
        # phi is efficieny of this reaction, number of molecules removed per event.
        rate[idx1:idx2] .= 4.0*pi*zeta*1.64e-4*(GRAIN_AREA_PER_H)*(1.0/mantle[dstep])*phi
    else
        rate[idx1:idx2] .= 0.0
    end

    # Desorption due to UV, partially from ISRF and partially from CR creating photons
    idx1=deuvcrReacs[1]
    idx2=deuvcrReacs[2]
    if desorb == 1 && uvdesorb == 1 && gama[j] <= ebmaxuvcr && mantle[dstep] >= 1.0e-20
        # 4.875d3 = photon flux, Checchi-Pestellini & Aiello (1992) via Roberts et al. (2007)
        # UVY is yield per photon.
        # 0.25 to change surface area of grain to cross-sectional area of grain (see freezeout)
        rate[idx1:idx2] = 0.25*GRAIN_AREA_PER_H*uv_yield*4.875e3*zeta*(1.0/mantle[dstep])
        # additional factor accounting for UV desorption from ISRF. UVCREFF is ratio of 
        # CR induced UV to ISRF UV.
        rate[idx1:idx2] = rate[idx1:idx2] .* (1+(radfield/uvcreff)*(1.0/zeta)*exp(-1.8*av[dstep]))
    else
        rate[idx1:idx2] = 0.0
    end

    idx1=thermReacs[1]
    idx2=thermReacs[2]
    if (thermdesorb == 1 && mantle[dstep] >= 1.0d-20) 
        for j in idx1:idx2
            # then try to overwrite with position in grain array
            for i in lbound(grainList,1):ubound(grainList,1)
                # See Cuppen, Walsh et al. 2017 review (section 4.1)
                if (grainList[i] ==  re1[j]) 
                    # Basic rate at which thermal desorption occurs
                    rate[j]=vdiff[i]*exp(-gama[j]/gasTemp[dstep])
                    # factor of 2.0 adjusts for fact only top two monolayers (Eq 8)
                    # becayse GRAIN_AREA_PER_H is per H nuclei, multiplying it by density gives area/cm-3
                    # that is roughly sigma_g.n_g from cuppen et al. 2017 but using surface instead of cross-sectional
                    # area seems more correct for this process.
                    rate[j]=rate[j]*(2.0/mantle[dstep])*SURFACE_SITE_DENSITY*GRAIN_AREA_PER_H*density[dstep]
                end
            end
        end
    else
        rate[idx1:idx2]=0.0
    end


    # Reactions on surface can be treated considering diffusion of reactants
    # See work of David Quenard 2017 Arxiv:1711.05184
    # abstracted to functions below for ease of reading
    idx1=diffReacs[1]
    idx2=diffReacs[2]
    for j in idx1:idx2
        rate[j]=diffusionReactionRate(j)
    end
    idx1=chemdesReacs[1]
    idx2=chemdesReacs[2]
    for j in idx1:idx2
        rate[j]=diffusionReactionRate(j)
    end

    # Basic gas phase reactions 
    # They only change if temperature has so we can save time with an if statement
    idx1=twobodyReacs[1]
    idx2=twobodyReacs[2]
    if lastTemp # = gasTemp[dstep] 
        rate[idx1:idx2] = alpha[idx1:idx2]*((gasTemp[dstep]/300.)^beta[idx1:idx2])*exp(-gama[idx1:idx2]/gasTemp[dstep]) 
    end

    lastTemp = gasTemp[dstep]
    if duplicates[1] # = 9999
        # this multiplies rate by 0 or 1 depending on whether gastemp>mintemp of a reaction
        rate[duplicates]=rate[duplicates]*min(real(floor(gasTemp[dstep]/minTemps)),1.0)
        # and this multiplies by 0,1 if gastemp>max temp
        rate[duplicates]=rate[duplicates]*min(real(floor(maxTemps/gasTemp[dstep])),1.0)
    end

    # Photoreactions for which we have a more detailed treatment
    h2dis = H2PhotoDissRate(h2Col,radField,av[dstep],turbVel) # H2 photodissociation
    rate[nrco] = COPhotoDissRate(h2Col,coCol,radField,av[dstep]) # CO photodissociation
    rate[nR_C_hv] = cIonizationRate(alpha[nR_C_hv],gama[nR_C_hv],gasTemp[dstep],ccol,h2col,av[dstep],radfield) # C photoionization
end

function diffusionReactionsRate(reacIndx)
    # want position of species in the grain array but gas phase species aren't in there
    # so store species index
    index1=re1[reacIndx]
    index2=re2[reacIndx]

    # then try to overwrite with position in grain array
    for i in lbound(grainList,1):ubound(grainList,1)
        if grainList[i] == index1 index1 = i end
        if grainList[i] == index2 index2 = i end
    end

    # Hasegawa 1992 diffusion rate. Rate that two species diffuse and meet on grain surface
    diffuseProb = vdiff[index1]*exp(-0.5*bindingEnergy[index1]/gasTemp[dstep])
    diffuseProb = diffuseProb+ (vdiff[index2]*exp(-0.5*bindingEnergy[index2]/gasTemp[dstep]))

    # probability a reactant will just desorb
    desorbProb = vdiff[index1]*exp(-bindingEnergy[index1]/gasTemp[dstep])
    desorbProb = desorbProb + vdiff[index2]*exp(-bindingEnergy[index2]/gasTemp[dstep]) 

    # Calculate classical activation energy barrier exponent
    reacProb = gama[reacIndx]/gasTemp[dstep]
    # Calculate quantum activation energy barrier exponent
    reducedMass = mass[grainList[index1]] * mass[grainList[index2]] / (mass[grainList[index1]] + mass[grainList[index2]])
    tunnelProb = 2.0 *CHEMICAL_BARRIER_THICKNESS/REDUCED_PLANCK * sqrt(2.0*AMU*reducedMass*K_BOLTZ*gama[reacIndx])

    # Choose fastest between classical and tunnelling
    if (reacProb.GT.tunnelProb) reacProb=tunnelProb end
    # set activationBarrier to probability of reaction Ruaud+2016
    reacProb=exp(-reacProb)

    # Overall reaction probability is chance of reaction occuring on meeting * diffusion rate
    reacProb = max(vdiff[index1],vdiff[index2]) * reacProb       


    #  Keff from Garrod & Pauly 2011 and Ruaud+2016
    #  Actual reaction probability is Preac/(Preac+Pevap+Pdiffuse), accounting for the other possible processes
    if DifFUSE_REACT_COMPETITION
       reacProb = reacProb/(reacProb + desorbProb + diffuseProb)
    end
    
    # see Eq A1 of Quenard et al. 2018
    # GAS_DUST_DENSITY_RATIO = 1/ frac abundance of dust so ratio/density is 1/number density of dust
    n_dust=density[dstep]/GAS_DUST_DENSITY_RATIO
    diffusionReactionRate=alpha[reacIndx] *reacProb* diffuseProb/(NUM_SITES_PER_GRAIN*n_dust)

    # Now adjust for fraction of this reaction's products that will desorb due to energy released
    if (reacIndx >= diffReacs[1] && reacIndx <= diffReacs[2]) 
        diffusionReactionRate = diffusionReactionRate * (1.0-desorptionFraction(reacIndx,index1,index2,true))
    else 
        diffusionReactionRate = diffusionReactionRate * desorptionFraction(reacIndx,index1,index2,false)
    end
end

function desorptionFraction(reacIndx,reactIndex1,reactIndex2,isDiff)
    # integer :: reacIndx,reactIndex1,reactIndex2,degreesOfFreedom
    # integer :: productIndex[4]
    # LOGICAL :: isDiff

    # double precision :: deltaEnthalpy,maxBindingEnergy,epsilonCd,productEnthalpy
    # double precision, parameter :: EFFECTIVE_SURFACE_MASS = 120.0

    
    # Get indices of grain surface version of products products 
    productIndex = 0
    productIndex = 0.0
    maxBindingEnergy = 0.0
    productEnthalpy = 0.0

    if isDiff
        for i in lbound(grainList,1):ubound(grainList,1)
            if (grainList[i]  == p1[reacIndx]) productIndex[1]=grainList[i]
            # Go through grain list and try to find species in product list
            # If it has a binding energy larger than largest energy found so far, update maxBindingEnergy
            if (grainList[i]  == p1[reacIndx]) 
                productIndex[1] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end

            if (grainList[i]  == p2[reacIndx]) 
                productIndex[2] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end

            if (grainList[i]  == p3[reacIndx]) 
                productIndex[3] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end

            if (grainList[i]  == p4[reacIndx]) 
                productIndex[4] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end

        end
    else
        for i in lbound(gasGrainList,1):ubound(gasGrainList,1)
            if (gasGrainList[i]  == p1[reacIndx]) 
                productIndex[1] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end
            if (gasGrainList[i]  == p2[reacIndx]) 
                productIndex[2] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end
            if (gasGrainList[i]  == p3[reacIndx]) 
                productIndex[3] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end
            if (gasGrainList[i]  == p4[reacIndx]) 
                productIndex[4] = grainList[i]
                productEnthalpy=productEnthalpy+formationEnthalpy[i]
                if (bindingEnergy[i] >= maxBindingEnergy) maxBindingEnergy=bindingEnergy[i]
            end
        end
    end
    
    # 
    # epsilonCd is the fraction of kinetic energy kept my the product when it collides with grain surface
    epsilonCd = mass(productIndex[1]) + mass(productIndex[2]) + mass(productIndex[3]) + mass(productIndex[4])
    epsilonCd = ((epsilonCd - EFFECTIVE_SURFACE_MASS) / (epsilonCd + EFFECTIVE_SURFACE_MASS))^2
    
    # Now calculate the change in enthalpy of the reaction.
    deltaEnthalpy= formationEnthalpy(reactIndex1)+formationEnthalpy(reactIndex2)-productEnthalpy
    
    # Convert from kcal to J, from J to K and from moles-1 to reactions-1
    deltaEnthalpy = deltaEnthalpy*4.184e03/(1.38054e-23*6.02214129e23)
    #  Total energy change includes activation energy of the reaction # 
    deltaEnthalpy = deltaEnthalpy + gama[reacIndx]


    if (deltaEnthalpy == 0.00) deltaEnthalpy = 1e-30 end

    # Degrees of freedom = 3 * number of atoms in the molecule
    degreesOfFreedom = atomCounts(productIndex[1])
    if (productIndex[2] != 0) degreesOfFreedom = max(degreesOfFreedom,atomCounts[productIndex[2]])
    if (productIndex[3] != 0) degreesOfFreedom = max(degreesOfFreedom,atomCounts[productIndex[3]])
    if (productIndex[4] != 0) degreesOfFreedom = max(degreesOfFreedom,atomCounts[productIndex[4]])                    
    degreesOfFreedom = 3 * degreesOfFreedom
        
    desorptionFraction = exp((-maxBindingEnergy*real(degreesOfFreedom)) / (epsilonCd * deltaEnthalpy))
    
    if deltaEnthalpy < 0       # < If reaction is endothermic, no CRD
        desorptionFraction = 0.
    end
    
    if GRAINS_HAVE_ICE
        desorptionFraction = desorptionFraction/10    # < See Minisalle et al. 2016 for icy grain surface.
        #  Special case of OH+H, O+H, N+N on ices, see same paper
        if (re1[reacIndx] == ngn &&re2[reacIndx] ==ngn) desorptionFraction = 0.5 end
        if ((re1[reacIndx] ==ngo &&re2[reacIndx] ==nh) || (re1[reacIndx] == nh &&re2[reacIndx] ==ngo)) desorptionFraction = 0.3 end
        if ((re1[reacIndx] ==ngoh &&re2[reacIndx] ==nh) || (re1[reacIndx] ==nh &&re2[reacIndx] ==ngoh)) desorptionFraction = 0.25 end
    end
end

end