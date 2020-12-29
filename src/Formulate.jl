using CSV 
using DataFrames

include("InitialNetworkConditions.jl")
include("ChemicalNetwork.jl")


function formulate(speciesFilepath::String, reactionsFilepath::String, icFilepath::String, p::Parameters, tspan::Tuple{Float64,Float64})
    reactionsData = CSV.read(reactionsFilepath, DataFrame)
    speciesData = CSV.read(speciesFilepath, DataFrame)
    initialConditionsData = CSV.read(icFilepath, DataFrame)

    calculateRates!(reactionsData, p)
    filterReactionData!(reactionsData, speciesData.name)

    inc = InitialNetworkConditions(icFilepath)
    fillInitialNetworkConditions(inc, speciesData.name)
    
    u0, network = createNetwork(inc, speciesData.name, reactionsData, p)
    ChemicalNetworkProblem(network, speciesData.name, u0, tspan)
end
