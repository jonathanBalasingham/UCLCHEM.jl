module Constants

const C  = 2.99792458D+10 #Speed of light in cgs
const K_BOLTZ = 1.38065040D-16 # Boltzmann constant cgs
const HP = 6.62606896D-27 #Planck constant in cgs
const REDUCED_PLANCK = 1.054571628d-27
const MH = 1.67262164D-24 #H nucleus mass in cgs
const AMU = 1.66053892d-24 #atomic mass unit in cgs
const K_BOLTZ_SI = 1.38d-23 #Boltzmann constant SI
const PC = 3.086d18 #parsec in cgs
const au = 2.063d5 #1 AU in cgs
const KM = 1.0d5 #kilometre in cgs
const SECONDS_PER_YEAR = 3.16d7
const T_CMB = 2.73
const EV = 1.60217646D-12 # electron volt in erg
const GRAV_G = 6.674d-8 #gravitational constant in cgs
const SB_CONST=5.6704d-5 #Stefan Boltzmann constant in cgs

function pair_insertion_sort!(array::T) where {T <: AbstractArray}
    for i in 2:2:length(array)-1
        t1 = min(array[i], array[i+1])
        t2 = max(array[i], array[i+1])
        j = i - 1
        while true
            array[j+2] = array[j]
            j =- 1
            j >= 1 && array[j] > t2 || break
        end
        array[j+2] = t2
        while true
           array[j+1] = array[j]
           j =- 1
           j >= 1 && array[j] > t1 || break
        end
        array[j+1] = t1
    end

    if last % 2 == 0
        t1 = array[last]
        for j in last-1:-1:1
            if array[j] <= t1 return end
            array[j+1] = array[j]
        end
        array(j+1)=t1
    end
end

export pair_insertion_sort!
export C, K_BOLTZ, K_BOLTZ_SI, KM, PC, au, SB_CONST, GRAV_G, EV, T_CMB, 
        SECONDS_PER_YEAR, AMU, HP, REDUCED_PLANCK, MH 


end