module Sputter
    function sputtering(abund::Array{Float64})
        # DOUBLE PRECISION :: abund(nspec+1,points)
        # INTENT(INOUT) :: abund
        # DOUBLE PRECISION :: sputterRate=0.0
        # INTEGER :: iSpec
  
        # loop over projectile species and get rates of change of mantle for each, summing them
        for iSpec in 1:size(projectiles) # # # #  Make projectiles array in initialize
              sputterRate += iceYieldRate(mass[projectiles[i],density[dstep]*abund[projectiles[i,dstep]]])
        end
  
        # Total rate is sputterRate (per grain) multiplied by grain number density
        
  
        # integrate that forward (check currentTime/targetTime) think I want to go
        # from currentTimeOld to currentTime.
  
        # Should get a portion of current mantle to remove, add that to gas phase.
        # eg abund(gasGrainList)+=sputterFrac*abund(grainList)
              # abund(grainList)-=sputterFrac*abund(grainList)
    end

    function iceYieldRate(projectileMass,projectileAbundance)
        # DOUBLE PRECISION projectileMass,projectileAbundance,sConst,driftVel
        iceBindingEnergy=0.53*1.6e-12
        # DOUBLE PRECISION :: lowerLimit,upperLimit
  
        sConst=(driftVel*driftVel)/(2.0*temp[dstep]*K_BOLTZ_CGS)
        sConst=sqrt(sConst)
  
        # eta is effectively reduced mass of the collision
        eta=4.0 * iceYieldEfficiency*projectileMass*targetMass*(projectileMass+targetMass)^(-2)
        epso=max(1.0,4.0 * eta)
  
        # Lower limit is xth in Jimenez-Serra et al. 2008
        lowerLimit=sqrt(epso*iceBindingEnergy/(eta*K_BOLTZ_CGS*temp[dstep]))
        # Upper limit is just where the integrand goes to zero
        upperLimit=iceYieldIntegralLimit(lowerLimit,projectileMass)
  
        iceYieldRate=trapezoidIntegrate(iceYieldIntegrand,lowerLimit,upperLimit)
        iceYieldRate=iceYieldRate*grainRadius*grainRadius*sqrt(8.0*K_BOLTZ_CGS*temp[dstep]*pi/projectileMass)
        iceYieldRate=iceYieldRate*projectileAbundance
    end

    function iceYieldIntegrand(x::T, projectileMass::T) where {T<:Float64}  
        yieldConst = 8.2e-4,
        iceYieldEfficiency=0.8 # 
        targetMass=18.0 # Mass of H2O in amu
        iceBindingEnergy=0.53*1.6e-12
  
        # this is s from exp(x+s), varies only with mass so should probably precalculate and save constant 
        s = sConst*sqrt(projectileMass)
        E = (x^2)*K_BOLTZ_CGS*temp
        eps = eta*E/iceBindingEnergy
        # this yield is for ice. There's a different one for cores (Appendix B Jimenez-Serra 2008)
        yield = yieldConst*(eps-epso)^2/(1 + (eps/30.)^(4/3))
        func_iceHe = yield*(x^2)*(exp(-(x-s)^2)-exp(-(x+s)^2))
        
    end

    function iceYieldIntegralLimit(xth::T, projectileMass::T) where {T<:Float64}
        i = 1
        while (iceYieldIntegrand(iceYieldIntegralLimit,projectileMass) < 1e-200)
              iceYieldIntegralLimit=xth+(1d3-xth)*(0.5^i)
        end   
        return iceYieldIntegralLimit     
    end
end