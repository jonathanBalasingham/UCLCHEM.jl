using Plots

function visualize(sol::ChemicalNetworkSolution; species=nothing, size=(1200,800))
    if isnothing(species)
        plot(sol.t ./ (3600 * 24 * 365), 
            log10.(transpose(reduce(hcat, sol.u))), 
            size=size, labels=permutedims(sol.species),
            xticks=([10^0, 10,10^2,10^3,10^4,10^5,10^6, 10^7, 10^8],["10^0", "10","10^2","10^3","10^4","10^5","10^6","10^7", "10^8"]),
            xaxis=:log10,
            legend=:bottomright)
        else
        inds =  indexin(species, sol.species)
        plot(sol.t ./ (3600 * 24 * 365), 
             log10.(transpose(reduce(hcat, sol.u))[:,inds]),
             size=size, labels=permutedims(sol.species[inds]),
             xticks=([10^0, 10,10^2,10^3,10^4,10^5,10^6, 10^7, 10^8],["10^0", "10","10^2","10^3","10^4","10^5","10^6","10^7", "10^8"]),
             xaxis=:log10,
             legend=:bottomright)
    end
end