using Flux
#include("../EchoStateNetwork.jl")


chain = Chain(Dense(4, 32), Dense(32, 32), Dense(32, 4))

function eval_model(x)
    out = chain(x)
    #Flux.reset!(chain)
    out
end

loss(x, y) = abs(sum((eval_model(x) .- y)))

ps = Flux.params(chain)

# use the ADAM optimizer. It's a pretty good one!
opt = Flux.ADAM()

data = (sol_as_matrix[:, 1:end-1]', sol_as_matrix[:, 2:end]')

println("Training loss before = ", sum(loss.(train_data, train_labels)))
#println("Test loss before = ", sum(loss.(test_data, test_labels)))

# callback function during training
num_epochs = 10
evalcb() = @show(sum(loss.(eachrow(data[1]), eachrow(data[2]))))
fake_data = Iterators.repeated((), 30)
loss1() = sum(loss.(eachrow(data[1]), eachrow(data[2])))


using Flux: @epochs
@epochs num_epochs Flux.train!((X) -> sum(loss(X[1], X[2])), ps, zip(eachrow(data[1]), eachrow(data[2])), opt, cb = Flux.throttle(evalcb, 1))