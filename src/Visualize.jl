using Plots

function visualize(sol::ChemicalNetworkSolution; species=nothing, size=(1300,800), take_every_nth=1)
    if isnothing(species)
        plot(sol.t[1:take_every_nth:end] ./ (3600 * 24 * 365), 
            log10.(transpose(reduce(hcat, sol.u[1:take_every_nth:end]))), 
            size=size, 
            labels=permutedims(sol.species),
            xticks=([10^0, 10,10^2,10^3,10^4,10^5,10^6],["10⁰", "10","10²","10³","10⁴","10⁵","10⁶"]),
            yticks=([10^-15, 10^-13, 10^-11, 10^-9, 10^-7, 10^-5, 10^-3, 10^-1, 10^1],["10⁻¹⁵", "10⁻¹³", "10⁻¹¹", "10⁻⁹", "10⁻⁷", "10⁻⁵", "10⁻³", "10⁻¹", "10¹"]),
            xaxis=:log10,
            legend=:outertopright)
        else
        inds =  indexin(species, sol.species)
        plot(sol.t[1:take_every_nth:end] ./ (3600 * 24 * 365), 
             log10.(transpose(reduce(hcat, sol.u[1:take_every_nth:end]))[:,inds]),
             #size=size, 
             labels=permutedims(sol.species[inds]),
             xticks=([10^0, 10,10^2,10^3,10^4,10^5,10^6],["10⁰", "10","10²","10³","10⁴","10⁵","10⁶"]),
             yticks=(log10.([10^-15, 10^-13, 10^-11, 10^-9, 10^-7, 10^-5, 10^-3, 10^-1, 10^1]),["10⁻¹⁵", "10⁻¹³", "10⁻¹¹", "10⁻⁹", "10⁻⁷", "10⁻⁵", "10⁻³", "10⁻¹", "10¹"]),
             xaxis=:log10,
             legend=:outertopright)
    end
    xlims!(1, 10^7)
    ylims!(log10.((10^-15, 10)))
    xlabel!("Time / Years")
    ylabel!("Xₛₚₑ")
end