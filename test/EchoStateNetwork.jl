using SparseArrays, LinearAlgebra
import Base.:*

struct EchoStateReservoir{F::Function, X<:AbstractArray, Y<:AbstractArray }
    W::X
    b::Y
    f::F
end

EchoStateReservoir(w::A, b::B, f::C) where {C<:Function, A<:AbstractArray, B<:AbstractArray} = EchoStateReservoir{C, A, B}(w,b,f)

function EchoStateReservoir(reservoir_size::Int, spectral_radius::Float64, sparsity::Float64; activation=tanh)
    reservoir_size < 1 ? error("Invalid reservoir size: $reservoir_size") : true
    0 < sparsity < 1 ? true : error("Sparsity out of range of (0,1)")
    r = sprand(reservoir_size, reservoir_size, sparsity) |> Matrix 
    EchoStateReservoir(r .*  spectral_radius / maximum(abs.(eigvals(r))), zeros(reservoir_size), activation)
end


function (res::EchoStateReservoir)(input_layer::AbstractArray)
    res.f(res * input_layer)
end

*(esr::EchoStateReservoir, input::AbstractArray) = esr.W * input


mutable struct EchoStateNetwork
    input_layer
    W
end

function (esn::EchoStateNetwork)(x::AbstractArray)

end

function train(esn::EchoStateNetwork, training_data::AbstractArray; method=:leastsquares)

    if method == :leastsquares

    else if method == :iterative

    else
        error("Invalid training method: $method")
    end
end
