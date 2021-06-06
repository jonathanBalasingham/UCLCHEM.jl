using CSV, DataFrames
using Flux
using Flux.Data

train = CSV.read("/home/jon/UCLCHEM/UCLCHEM.jl/test/PoC/data/roberson_weights.csv", DataFrame, header=false) |> Matrix
X = train[1:2:7000, 2:4]
y = train[1:2:7000, 5:end]

#train_loader = DataLoader((X', y), batchsize=10, shuffle=true)
println("Making chain")
chain = Flux.Chain(Dense(size(X, 2),256, tanh), Dense(256, 512, tanh), Dense(512, 1024, tanh),  Dense(1024, size(y, 2), tanh))

function predict_err(xk, yk)
    pred = chain(xk)
    Flux.Losses.mae(pred, yk)
end

#= Full batch =#
loss() = sum(predict_err.(eachrow(X),eachrow(y)))
#=
function loss()
    sample_index = rand(1:size(X, 1)-200)
    x_subsample = X[sample_index:sample_index+200, :]
    y_subsample = y[sample_index:sample_index+200, :]
    sum(predict_err.(eachrow(x_subsample),eachrow(y_subsample)))
end
=#
parameters = params(chain)
opt = ADAM(0.001)
data = Iterators.repeated((), 50)

evalcb() = @show(loss())
epochs = 2

for i in 1:epochs
    print("Starting loss at epoch $i is: ")
    evalcb()
    @time Flux.train!(loss, parameters, data, opt)
    print("Ending loss at epoch $i is: ")
    evalcb()
end

using BSON: @save, @load
@save "roberson_network.bson" chain
