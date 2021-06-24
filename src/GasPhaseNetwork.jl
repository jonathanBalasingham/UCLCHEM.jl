using Catalyst, DataFrames, CSV


include("Reaction.jl")
include("Rates.jl")
include("InitialNetworkConditions.jl")
include("ChemicalNetwork.jl")


function read_in_reactions(reaction_file::String)
    d = CSV.read(reaction_file, DataFrame)
    d.re3 = d.re3 .|> string
    d.prod4 = d.prod4 .|> string
    d
end

function get_species_from_reactions(reaction_data::DataFrame)
    species = []
    for col in ["re1", "re2", "re3", "prod1", "prod2", "prod3", "prod4"]
        species = unique(skipmissing([species; reaction_data[!, col]]))
    end
    remove_these = ["NAN", "CRP", "PHOTON", "CRPHOT", "NaN", "nan", NaN, 
                    "DIFF", "FREEZE", "THERM", "CHEMDES", "DESCR", "DESOH2", "DEUVCR", "H2FORM"]
    filter!(x -> !(x in remove_these), species)
    filter!(!isnan, species)
    filter!(!ismissing, species)
    return Set(species)
end

function create_gas_phase_network(reactions::DataFrame; species=nothing)

    if isnothing(species)
        species = get_species_from_reactions(reactions) .|> Symbol
    else
        species = species .|> Symbol
    end

    network = make_empty_network()
    rates = reactions[!, end]
    #@parameters t rates[1: length(rates)]
    addparam!(network, (@parameters t)[1])
    make_variable = x -> (@variables $x(t))[1]

    for spec in species
       addspecies!(network, (@variables $spec(t))[1]) 
    end
    remove_these = ["NAN", "CRP", "PHOTON", "CRPHOT", "NaN", "nan", NaN, 
                    "DIFF", "FREEZE", "THERM", "CHEMDES", "DESCR", "DESOH2", "DEUVCR", "H2FORM"]

    for reaction_row in eachrow(reactions)
        re = filter(x -> !ismissing(x) && !isnan(x) && !(x in remove_these), reaction_row[1:3] |> Array) .|> Symbol
        pr = filter(x -> !ismissing(x) && !isnan(x) && !(x in remove_these), reaction_row[4:7] |> Array) .|> Symbol
        rate = reaction_row[end]
        rx = Reaction(rate, make_variable.(re), make_variable.(pr))
        addreaction!(network, rx)
    end
    network
end


function formulate_ode_system(rx_network)
    sys = convert(ODESystem, rx_network)
    ODESystem(sys.eqs, (@parameters t)[1])
end


function formulate_ode_problem(rx_network, u0, tspan)
    sys = formulate_ode_system(rx_network)
    u0map = map((x,y) -> Pair(x,y), species(rx_network), u0)
    ODEProblem(rx_network, u0map, tspan)
end


function read_in_initial_conditions(icfp::String)
    ic = CSV.read(icfp, DataFrame)
    return Dict(zip(ic.species, ic.concentration))
end


function get_species_names(network)
    species(network) .|> 
        string .|>
        x -> replace(x, "(t)" => "") .|>
        x -> replace(x, "\"" => "") .|>
        x -> replace(x, "var" => "")
end


function create_u0(network::ReactionSystem, initial_conditions::Dict; default_value=10e-20)
    species_names = species(network) .|> 
        string .|>
        x -> replace(x, "(t)" => "") .|>
        x -> replace(x, "\"" => "") .|>
        x -> replace(x, "var" => "")
    
    u0 = repeat([default_value], length(species_names))
    for (i, sn) in enumerate(species_names)
        if sn in keys(initial_conditions)
            u0[i] = initial_conditions[sn]
        end
    end
    u0
end

function formulate_all(rfp::String, icfp::String, p; tspan=(0., 3600. *24. *365. *1000000))
    ics = read_in_initial_conditions(icfp)
    reactions_data = read_in_reactions(rfp)
    calculateRates!(reactions_data, p)
    network = create_gas_phase_network(reactions_data)
    sys = formulate_ode_system(network)
    species_name = get_species_names(network)
    u0 = create_u0(network, ics)
    cnp = ChemicalNetworkProblem(sys, species_name, u0, tspan)
    return cnp
end
