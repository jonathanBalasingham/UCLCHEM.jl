<<<<<<< HEAD
using SparseArrays, LinearAlgebra
import Base.:*
using Flux

struct EchoStateReservoir{F, X<:AbstractArray, Y<:AbstractArray}
    W::X
    b::Y
    f::F
end

EchoStateReservoir(w::A, b::B, f::C) where {C<:Function, A<:AbstractArray, B<:AbstractArray} = EchoStateReservoir{C, A, B}(w,b,f)

function EchoStateReservoir(reservoir_size::Int, spectral_radius::Float64, sparsity::Float64; activation=identity)
    reservoir_size < 1 ? error("Invalid reservoir size: $reservoir_size") : true
    0 < sparsity < 1 ? true : error("Sparsity out of range of (0,1)")
    
    r = sprand(reservoir_size, reservoir_size, sparsity) |> Matrix 
    EchoStateReservoir(r .*  spectral_radius / maximum(abs.(eigvals(r))), zeros(reservoir_size), activation)
end


function (res::EchoStateReservoir)(input_layer::AbstractArray)
    res.f(res * input_layer)
end

*(esr::EchoStateReservoir, input::AbstractArray) = esr.W * input


#=
This could be changed to a Flux Chain, but I need to figure out
a way to make the reservoir to behave recurrently instead of as 
just a normal feed-forward layer.
=#
mutable struct EchoStateNetwork
    input_layer::Dense
    W::EchoStateReservoir
    output_layer::Dense
    state::AbstractArray
end

function EchoStateNetwork(input_layer::Dense, W::EchoStateReservoir, output_layer::Dense)
    input_size = size(input_layer.b, 1)
    res_size = size(W.b, 1)
    viable_res_size = Int(floor(approx_res_size/input_size)*input_size)
    if viable_res_size != res_size
        error("Input Layer and reservoir do not have compatible size, reservoir needs size: $viable_res_size")
    end

    new(input_layer, W, output_layer, zeros(size(W.b, 1)))
end


function (esn::EchoStateNetwork)(x::AbstractArray)
    input_layer_output = x |> esn.input_layer
    output_layer * esn.W.f(esn.W * esn.state + input_layer_output)
end

function train(esn::EchoStateNetwork, train_sequence::AbstractArray; method=:leastsquares, beta=0.01)

    if method == :leastsquares
        res_output = esn(train_sequence)
        term2 = (res_output*res_output') .+ beta
        esn.output_layer = (train_sequence*res_output') * inv(term2)
    elseif method == :iterative

    else
        error("Invalid training method: $method")
    end
end
||||||| merged common ancestors
=======
using SparseArrays, LinearAlgebra
import Base.:*
using Flux

struct EchoStateReservoir{F, X<:AbstractArray, Y<:AbstractArray}
    W::X
    b::Y
    f::F
end

EchoStateReservoir(w::A, b::B, f::C) where {C<:Function, A<:AbstractArray, B<:AbstractArray} = EchoStateReservoir{C, A, B}(w,b,f)

function EchoStateReservoir(reservoir_size::Int, spectral_radius::Float64, sparsity::Float64; activation=identity)
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
    input_layer::Dense
    W::EchoStateReservoir
    output_layer::Dense
    state::AbstractArray
end

function EchoStateNetwork(input_layer::Dense, W::EchoStateReservoir, output_layer::Dense)
    input_size = size(input_layer.b, 1)
    res_size = size(W.b, 1)
    viable_res_size = Int(floor(approx_res_size/input_size)*input_size)
    if viable_res_size != res_size
        error("Input Layer and reservoir do not have compatible size, reservoir needs size: $viable_res_size")
    end

    new(input_layer, W, output_layer, zeros(size(W.b, 1)))
end


function (esn::EchoStateNetwork)(x::AbstractArray)
    input_layer_output = x |> esn.input_layer
    output_layer * esn.W.f(esn.W * esn.state + input_layer_output)
end

function train(esn::EchoStateNetwork, training_data::AbstractArray; method=:leastsquares)

    if method == :leastsquares

    elseif method == :iterative

    else
        error("Invalid training method: $method")
    end
end
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77
