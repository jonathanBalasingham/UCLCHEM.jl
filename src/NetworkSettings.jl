include("Rates.jl")

struct NetworkSettings
    speciesFilepath::String
    reactionsFilepath::String
    initialConditionsFilepath::String
    p::Parameters
end