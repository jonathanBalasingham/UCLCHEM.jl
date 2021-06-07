module PhotoReactions

UV_FAC = 3.02
#Below are arrays for self-shielding of CO and H2
start = true
NUM_LAMBDA = 30

LAMBDA_GRID = [910.0 950.0 1000.0 1050.0 1110.0
               1180.0 1250.0 1390.0 1490.0 1600.0 
               1700.0 1800.0 1900.0 2000.0 2100.0 
               2190.0 2300.0 2400.0 2500.0 2740.0 
               3440.0 4000.0 4400.0 5500.0 7000.0 
               9000.0 12500.0 22000.0 34000.0 1.0e9]

XLAMBDA_GRID = [5.76 5.18 4.65 4.16 3.73  
                3.40 3.11 2.74 2.63 2.62  
                2.54 2.50 2.58 2.78 3.01  
                3.12 2.86 2.58 2.35 2.00  
                1.58 1.42 1.32 1.00 0.75  
                0.48 0.28 0.12 0.05 0.00]

XLAMBDA_DERIV = zeros(30)

startr = true

#  12CO line shielding data from van Dishoeck  Black (1988, ApJ, 334, 771, Table 5)
DIMCO = 7 
DIMH2 = 6

NCO_GRID = [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0]
NH2_GRID = [18.0, 19.0, 20.0, 21.0, 22.0, 23.0]

SCO_GRID = [0.000  -1.408e-02 -1.099e-01 -4.400e-01 -1.154  -1.888  -2.760  -4.001   
            -8.539e-02 -1.015e-01 -2.104e-01 -5.608e-01 -1.272  -1.973  -2.818  -4.055   
            -1.451e-01 -1.612e-01 -2.708e-01 -6.273e-01 -1.355  -2.057  -2.902  -4.122   
            -4.559e-01 -4.666e-01 -5.432e-01 -8.665e-01 -1.602  -2.303  -3.146  -4.421   
            -1.303  -1.312  -1.367  -1.676  -2.305  -3.034  -3.758  -5.077   
            -3.883  -3.888  -3.936  -4.197  -4.739  -5.165  -5.441  -6.446]
SCO_DERIV = zeros(8, 6)

function H2PhotoDissRate(NH2, radField, av, turbVel)
    baseRate = 5.18e-11
    xl = 1000.0 
    radWidth=8.0e7 

    dopplerWidth = turbVel/(xl*1.0d-8)
    baseRate * (radField/1.7) * scatter(xl, av) * H2SelfShielding(NH2, dopplerWidth, radWidth)
end

function COPhotoDissRate(NH2,NCO,radField,av)
    ssf = COSelfShielding(NH2, NCO)
    lba = lbar(NCO, NH2)
    sca = scatter(lba, av)

    # The reason why rad is divided by 1.7 is that the alphas are for Draine and the rad is in 
    # Habing units
    COPhotoDissRate = (2.0e-10) * (radfield/1.7) * ssf * sca
end

function cIonizationRate(alpha,gamma,gasTemp,NC,NH2,av,radfield)
    #  Calculate the optical depth in the CI absorption band, accounting
    #  for grain extinction and shielding by CI and overlapping H2 lines
    TAUC = gamma*av+1.1e-17*NC+(0.9*gasTemp^0.27*(NH2/1.59e21)^0.45)
    #  Calculate the CI photoionization rate
    RATE = alpha*(radfield/1.7)*exp(-TAUC)
end

function H2SelfShielding(NH2,dopplerWidth,radWidth)
    FPARA = 0.5
    FOSC = 1.0e-2
    taud  = FPARA * NH2 * 1.5e-2 * FOSC / dopplerWidth
     
    # calculate doppler contribution of self shielding function sj
    if (taud .eq. 0.0) 
       sj = 1.0
    elseif (taud .lt. 2.0) 
       sj = exp(-0.6666667*taud)
    elseif (taud .lt. 10.0) 
       sj = 0.638*taud^(-1.25)
    elseif (taud .lt. 100.0) 
       sj = 0.505*taud^(-1.15)
    else
       sj = 0.344*taud^(-1.0667)
    end

    # calculate wing contribution of self shielding function sr
    # if (taud.lt.0.0) taud=0.0 end
    if (radWidth .eq. 0.0)
       sr = 0.0
    else
       r  = radWidth/(1.7724539*dopplerWidth)
       t  = 3.02 * ((r*1.0d+03)^(-0.064))
       u  = sqrt(taud*r)/t
       sr = r/(t*sqrt(0.78539816+u^2.0))
    end
    
    # calculate total self shielding function fgk
    H2SelfShielding = sj + sr
end

function scatter(x1, av)
    c = [1.0, 2.006, -1.438, 7.364e-01, -5.076e-01, -5.920e-02]
    k1 = [7.514e-01, 8.490e-01, 1.013, 1.282, 2.005, 5.832]

    # optical depth in the visual
    tv = av/ 1.086
      
    # make correction for the wavelength considered
    tl = tv * xlambda(x1)
       
    # calculate attuentuation  due to dust scattering
    scatter = 0.0
    if tl < 1.0
        expo = k1[1] * tl
        if expo < 100.0 
            scatter = c[1] * exp(-expo)
        end
    else
        for i in 2:6
            expo = k1[i]*tl
            if expo < 100.0 
                scatter += c[i]*exp(-expo)
            end
        end
    end
    return scatter
end

function xlambda(LAMBDA)

    if START
        SPLINE(LAMBDA_GRID, XLAMBDA_GRID, NUM_LAMBDA,1.0e30,1.0e30, XLAMBDA_DERIV)
    end

    LAMBDA_VALUE = LAMBDA
    if LAMBDA < LAMBDA_GRID[1] LAMBDA_VALUE = LAMBDA_GRID[1] end
    if LAMBDA > LAMBDA_GRID[NUM_LAMBDA] LAMBDA_VALUE = LAMBDA_GRID[NUM_LAMBDA] end

    SPLINT(LAMBDA_GRID,XLAMBDA_GRID,XLAMBDA_DERIV,NUM_LAMBDA,LAMBDA_VALUE,xlambda)
    if XLAMBDA < 0.0 XLAMBDA = 0.0 end

end

function COSelfShielding(NH2,NCO)
    if startr
        splie2(NCO_GRID,NH2_GRID,SCO_GRID,DIMCO,DIMH2,SCO_DERIV)
        startr = false
    end
   
    lognco = log10(NCO+1.0)
    lognh2 = log10(NH2+1.0)
       
    if lognco < NCO_GRID[1]      lognco = NCO_GRID[1] end
    if lognh2 < NH2_GRID[1]      lognh2 = NH2_GRID[1] end
    if lognco > NCO_GRID[DIMCO]  lognco = NCO_GRID[DIMCO] end
    if lognh2 > NH2_GRID[DIMH2]  lognh2 = NH2_GRID[DIMH2] end
       
    splin2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,DIMCO,DIMH2,lognco,lognh2,COSelfShielding)
    COSelfShielding = 10.0^COSelfShielding
end

function lbar(u,w)
    # calculate lambda bar (in a) according to equ. 4 of van dishoeck
    # and black, apj 334, p771 (1988)
    #  --------------------------------------------------------------
    #        i/o parameter
    #        u : co column density in (cm-2)
    #        w : h2 column density in (cm-2)
    
    #         program variables
    #         lu : log10(co column density in cm-2)
    #         lw : log10(h2 column density in cm-2)
    
    # --------------------------------------------------------------
    lu = log10(abs(u)+1.0)
    lw = log10(abs(w)+1.0)
    
    lbar = (5675.0 - 200.6*lw) - (571.6 - 24.09*lw)*lu + (18.22 - 0.7664*lw)*lu^2
       
    # lbar represents the mean of the wavelengths of the 33
    # dissociating bands weighted by their fractional contribution
    # to the total rate of each depth. lbar cannot be larger than
    # the wavelength of band 33 (1076.1a) and not be smaller than
    # the wavelength of band 1 (913.6a).
    if lbar > 1076.1  lbar = 1076.1 end
    if lbar < 913.6  lbar =  913.6 end
end

function splie2(x1a,x2a,ya,m,n,y2a)
     # given an m by n tabulated function ya, and tabulated indepen-
     # dent variables x1a (m values) and x2a (n values), this routine
     # constructs one-dimensional natural cubic splines of the rows
     # of ya and returns the second-derivatives in the array y2a.
     # (copied from numerical recipes)
    
     # --------------------------------------------------------------
     # i/o parameter and program variables
     #         double precision  x1a(m), x2a(n), ya(m,n), y2a(m,n), ytmp(nn), y2tmp(nn)
     # --------------------------------------------------------------
     nn = 100
     y2tmp = zeros(nn)
     for i in 1:m
        for k in 1:n
            ytmp[k] = ya[j, k]
        end
        spline(x2a,ytmp,n,1.0e30,1.0e30,y2tmp)
        for k in 1:n
            y2a[j, k] = y2tmp[k]
        end
     end
end

function spline(x,y,n,yp1,ypn,y2)
    if yp1 >= 1.0e30
     # the lower boundary condition is set either to be
     # "natural"
        y2[1] =  0.0
        u[1]  =  0.0
    else
     # or else to have a specified first derivative.
        y2[1] = -0.5
        u[1]  = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
    end
   
     # this is the decomposition loop of the tridiagonal algorithm.
     # y2 and u are used for temporary storage of decomposed factors.
    for  i in 2:n-1
        sig   = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p     = sig*y2[i-1] + 2.0
        y2[i] = (sig-1.0)/p
        u[i]  = (6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p
    end
       
    if ypn >= 1.0e30
     #      the upper boundary condition is set either to be
     #      "natural"
        qn = 0.0
        un = 0.0
    else
     #      or else to have a specified first derivative.
        qn = 0.5
        un = (3.0/(x[n]-x[n-1])) *(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]))
    end
       
    y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0)
       
     # this is the backsubstitution loop of the tridiagonal algorithm
    for k in n-1:-1:1
        y2[k] = y2[k]*y2[k+1]+u[k]
    end
end


end