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

end