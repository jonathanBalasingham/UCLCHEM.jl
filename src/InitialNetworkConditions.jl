
struct InitialNetworkConditions
    initialConcentrations::Dict{String, Float64}
    settingsFilepath::String
end

function InitialNetworkConditions(fp::String)
    ics = CSV.read(fp, DataFrame)
    d = Dict{String, Float64}()
    map(row -> d[strip(row.species)] = row.concentration, eachrow(ics))
    InitialNetworkConditions(d, fp)
end

function fillInitialNetworkConditions(inc::InitialNetworkConditions, speciesList)
    for species in speciesList
        if !(species in keys(inc.initialConcentrations))
            inc.initialConcentrations[species] = 10^-20
        end
    end
end
