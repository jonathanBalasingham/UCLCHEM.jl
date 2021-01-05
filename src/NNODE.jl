using CUDA
using Flux, Zygote, DiffEqSensitivity, ForwardDiff, Random, Distributions
using DiffEqFlux, Adapt
using ModelingToolkit
using GalacticOptim
using Quadrature
using RuntimeGeneratedFunctions
import Tracker, Optim
import ModelingToolkit: value, nameof, toexpr, build_expr

abstract type NeuralPDEAlgorithm <: DiffEqBase.AbstractODEAlgorithm end

struct NNODE{C,O,P,K} <: NeuralPDEAlgorithm
    chain::C
    opt::O
    initθ::P
    autodiff::Bool
    kwargs::K
end

function NNODE(chain,opt=Optim.BFGS(),init_params = nothing;autodiff=false,kwargs...)
    if init_params === nothing
        if chain isa FastChain
            initθ = DiffEqFlux.initial_params(chain)
        else
            initθ,re  = Flux.destructure(chain)
        end
    else
        initθ = init_params
    end
    NNODE(chain,opt,initθ,autodiff,kwargs)
end