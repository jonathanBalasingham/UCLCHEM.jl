using Plots




function visualize(sol::ChemicalNetworkSolution; species=nothing, size=(1200,800))
    if isnothing(species)
        plot(sol.t ./ (3600 * 24 * 365), log10.(transpose(reduce(hcat, sol.u) .+ 1.)), size=size, labels=permutedims(sol.species))
    else
        inds =  indexin(species, sol.species)
        plot(sol.t ./ (3600 * 24 * 365), log10.(transpose(reduce(hcat, sol.u) .+ 1.))[:,inds], size=size, labels=permutedims(sol.species[inds]))
    end
end