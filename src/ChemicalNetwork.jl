include("Reaction.jl")
using ModelingToolkit


function prod(ind, species, reactions, Y)
    name = species[ind]
    # find reactions that have name in prod 1,2,3 or 4
    rows = filter(row -> name in row[4:7], reactions)
    p = 0
    for r in eachrow(rows)
        re1,re2,re3 = r[1:3]
        term = 1
        for (i,s) in enumerate(species)
            if s == re1
                term *= Y[i]
            elseif s == re2
                term *= Y[i]
            elseif s == re3
                term *= Y[i]
            end
        end
        p += term*r.rate
    end
    p
end

function loss(ind, species, reactions, Y)
    name = species[ind]
    # find reactions that have name in prod 1,2,3 or 4
    rows = filter(row -> name in row[1:3], reactions)
    # filter the row
    p = 0
    for r in eachrow(rows)
        pr1,pr2,pr3,pr4 = r[4:7]
        term = -1
        for (i,s) in enumerate(species)
            if s == pr1
                term *= Y[i]
            elseif s == pr2
                term *= Y[i]
            elseif s == pr3
                term *= Y[i]
            elseif s == pr4
                term *= Y[i]
            end
        end
        p += term*r.rate
    end
    p
end

function getHydrogenIndex(slist)
    i = 1
    for sp in slist
        if sp == "H"
            return i 
        else
            i += 1
        end
    end
    println("Couldn't find H, not adding H2 formation")
end

function createNetwork(IC::InitialNetworkConditions, species::Array{String}, reactionsData::DataFrame, p::Parameters)
    ModelingToolkit.@parameters t
    @variables y[1:length(species)](t)
    @derivatives D'~t

    nH = 1. # CHANGE ME
    H_index = getHydrogenIndex(species)
    h2_rate = 10^-17*sqrt(p.T)*y[H_index]
    ydot = D.(y)
    eqs = []
    u0 = Float64[]
    for (i,yt) in enumerate(y)
        name = species[i]
        push!(u0, IC.initialConcentrations[name])
        if name == "H"
            push!(eqs, ydot[i] ~ prod(i,species,reactionsData, y) + y[i]*loss(i,species,reactionsData,y) - h2_rate)
        elseif name == "H2"
            push!(eqs, ydot[i] ~ prod(i,species,reactionsData, y) + y[i]*loss(i,species,reactionsData,y) + h2_rate)
        else
            push!(eqs, ydot[i] ~ prod(i,species,reactionsData, y) + y[i]*loss(i,species,reactionsData,y))
        end
    end
    eqs = eqs .|> simplify
    u0, ODESystem(eqs, name=:uclchem)
end

struct ChemicalNetworkProblem
    network::ODESystem
    species::Array{String,1}
    u0::Array{Float64,1}
    tspan::Tuple{Float64, Float64}
end

struct ChemicalNetworkSolution
    t
    u::Array{Array{Float64,1},1}
    species::Array{String,1}
end

species(network::ChemicalNetworkProblem) = network.species