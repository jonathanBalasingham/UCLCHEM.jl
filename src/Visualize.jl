using Plots

function visualize(sol::ODESolution, speciesList, size=(1200,800))
    plot(sol, log10.(sol.u[1] .+ 1.), size=size, labels=permutedims(speciesList))
end