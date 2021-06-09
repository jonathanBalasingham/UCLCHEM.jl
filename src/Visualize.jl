using Plots

function visualize(sol::ChemicalNetworkSolution; species=nothing, size=(1200,800), take_every_nth=10)
    if isnothing(species)
        plot(sol.t[1:take_every_nth:end] ./ (3600 * 24 * 365), 
            log10.(transpose(reduce(hcat, sol.u[1:take_every_nth:end]))), 
            size=size, labels=permutedims(sol.species),
            xticks=([10^0, 10,10^2,10^3,10^4,10^5,10^6, 10^7, 10^8],["10^0", "10","10^2","10^3","10^4","10^5","10^6","10^7", "10^8"]),
            xaxis=:log10,
            legend=:outertopright)
        else
        inds =  indexin(species, sol.species)
        plot(sol.t[1:take_every_nth:end] ./ (3600 * 24 * 365), 
             log10.(transpose(reduce(hcat, sol.u[1:take_every_nth:end]))[:,inds]),
             size=size, labels=permutedims(sol.species[inds]),
             xticks=([10^0, 10,10^2,10^3,10^4,10^5,10^6, 10^7, 10^8],["10^0", "10","10^2","10^3","10^4","10^5","10^6","10^7", "10^8"]),
             xaxis=:log10,
             legend=:outertopright)
    end
end