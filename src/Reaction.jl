import Base.isnan

include("./Rates.jl")

isnan(s::String) = s == "NaN" || s == "NAN"
isspecies(s::String) = !(s in ["PHOTON", "CRP", "CRPHOT"])

function reactionCriteria(reaction, species)
    reAndProd = filter(!ismissing, Array(reaction[1:7]))
    reAndProd = filter(!isnan, reAndProd)
    reAndProd = filter(isspecies, reAndProd)
    for i in reAndProd
        if !(i in species)
            return false
        end
    end
    true
end

function filterReactionData!(rData, species)
    filter!(row -> reactionCriteria(row, species), rData)
end


function createReactionList(reactionData, includedSpeciesNames)
    reactionList = []

    for row in eachrow(reactionData)

        reactantsAndProducts = filter(x -> !isnan(x), Array(row[1:7]))
        include = true

        for sp in reactantsAndProducts
            if !(sp in includedSpeciesNames)
                include = false
            end
        end

        if include
            re = filter(x -> !isnan(x), [row.re1, row.re2, row.re3])
            prod = filter(x -> !isnan(x), [row.prod1, row.prod2, row.prod3, row.prod4])

            push!(reactionList, Reaction(re, prod, row.alpha, row.beta, row.gamma, row.tmin, row.tmax))
        end
    end
    reactionList
end