
R_AB(alpha::Float64, beta::Float64, gamma::Float64, T::Float64) = alpha*(T/(300))^beta*exp(-gamma/T)
R_CRProton(alpha::Float64, zeta::Float64) = alpha*zeta
R_CRPhoton(alpha::Float64, beta::Float64, gamma::Float64, zeta::Float64, omega::Float64, T::Float64, E::Float64) =
     alpha*gamma*(T/(300)^beta *((E*zeta)/(1-omega)))
R_UV(alpha::Float64, F_UV::Float64, gamma::Float64, A_v::Float64) = alpha*F_UV*exp(-gamma*A_v)
H2_form(T::Float64) = 10^-17*sqrt(T)

struct Parameters
    zeta::Float64
    omega::Float64
    T::Float64
    F_UV::Float64
    A_v::Float64
    E::Float64
    density::Float64
end

# per day
function calculateRates!(rdata, parameters)
    h2_formation_row = DataFrame(re1="H", re2="H", re3="H2FORM", 
                                 prod1="H2", prod2="NaN", prod3="NaN", prod4="NaN",
                                 alpha=0.0, beta=0.0, gamma=0.0, tmin=0.0, tmax=0.0)
    append!(rdata, h2_formation_row)
    rdata[!, "rate"] .= 0.0
    for row in eachrow(rdata)
        if row.re2 == "CRP"
            row.rate = R_CRProton(row.alpha,parameters.zeta)
        elseif row.re2 == "CRPHOT"
            # E parameter is missing!!
            row.rate = R_CRPhoton(row.alpha, row.beta, row.gamma, parameters.zeta, parameters.omega, parameters.T,1.0)
        elseif row.re2 == "PHOTON"
            row.rate = R_UV(row.alpha, parameters.F_UV, row.gamma, parameters.A_v)
        elseif row.re3 == "H2FORM"
            row.rate = H2_form(parameters.T)
        else
            row.rate = R_AB(row.alpha,row.beta,row.gamma,parameters.T) * parameters.density
        end
    end
end