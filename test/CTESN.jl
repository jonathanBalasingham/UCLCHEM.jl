#using Flux
using UCLCHEM
using Sundials
using Dates
using OrdinaryDiffEq
using ReservoirComputing
using ParameterizedFunctions
using Plots
using SQLite
using DelimitedFiles


function read_in_train_data(filepath::AbstractString, read_every::Int, include_time=false)
    data = Vector{Vector{Float64}}(undef, 0)
    file = open(filepath)
    counter = 0
    for i in eachline(file)
        if (counter % read_every == 0)
            b = IOBuffer(i)
            chunk = readdlm(b, ',', Float64)
            if length(data) == 0
                data = chunk
            else
                data = vcat(data, chunk)
            end
        end
        counter += 1
    end
    close(file)
    println(counter)
    if !include_time
        Matrix(data[:, 2:end]')
    else
        Matrix(data')
    end
end

function read_first_rows(filepath::AbstractString, read_first::Int, include_time=false)
    data = Vector{Vector{Float64}}(undef, 0)
    file = open(filepath)
    counter = 0
    for i in eachline(file)
        if (counter <= read_first)
            b = IOBuffer(i)
            chunk = readdlm(b, ',', Float64)
            if length(data) == 0
                data = chunk
            else
                data = vcat(data, chunk)
            end
        else
            break
        end
        counter += 1
    end
    close(file)
    println(counter)
    if !include_time
        Matrix(data[:, 2:end]')
    else
        Matrix(data')
    end
end


function create_and_test_esn(res_size,
                            radius,
                            degree,
                            activation,
                            sigma,
                            beta,
                            alpha, # leaking factor
                            nla_type,
                            extended_states,
                            train,
                            test,
                            predict_len;
                            make_plot=false)

    esn = ESN(  res_size,
                train,
                degree,
                radius,
                activation = activation, #default = tanh
                alpha = alpha, #default = 1.0
                sigma = sigma, #default = 0.1
                nla_type = nla_type, #default = NLADefault()
                extended_states = extended_states #default = false
    )

    @time W_out = ESNtrain(esn, beta)
    @time output = ESNpredict(esn, predict_len, W_out)
     
    test_error = sum(abs.(output - test))

    if test_error > 100
        return
    end

    db = SQLite.DB("./test/ctesn_params.sqlite")

    query = """
        INSERT INTO results VALUES (
            $res_size ,
            $radius ,
            $degree ,
            $sigma ,
            $beta , 
            $alpha ,
            $test_error
        )
    """

    DBInterface.execute(db, query) 
    DBInterface.close!(db)
    #print(test_error)
    file_name = "ESN_$res_size" * "_$radius" * "_$degree" * "_$sigma" * "_$beta" * "_$alpha" * "_.png"

    if make_plot
        p = plot(transpose(output[1:12,:]),layout=(4,3), label="predicted",size=(1200,800))
        plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
        savefig(p, "./test/CTESN_Plots/" * file_name)
    end
    print(file_name)
end


function main()
    res_sizes = [1001]
    radii = collect(1.1:.1:10)
    degrees = [10, 50]
    sigmas = collect(2.0:.2:20.0)
    betas = collect(0.0:.2:1.0)
    alphas = collect(.99:-.1:.49)


    rfp = "./test/input/reactions_postNN.csv"
    sfp = "./test/input/species_postNN.csv"
    icfp = "./test/input/initcond1.csv"

<<<<<<< HEAD
    T=10. 
    zeta = 1.
    omega = 0.5
    F_UV=1.
    A_v=2.
    E = 0.5
    dens = 1e4
    p = UCLCHEM.Parameters(zeta, omega, T, F_UV, A_v, E, dens)

    #tspan = (0., 10^7 * 365. * 24. * 3600.)
    tspan = (0., 10^6 * 365. * 24. * 3600.)

    nw_prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan)
    #prob = ODEProblem(nw_prob.network, nw_prob.u0, tspan)
    #sol = solve(prob, CVODE_BDF(), saveat = 24*360000*365)
    sol2 = UCLCHEM.solve(nw_prob, time_factor_pre_1000_years=20)
    println("Problem solved with CVODE")

    v = sol2.u
    data = Matrix(hcat(v...))
    shift = 1
    train_len = 100000
    predict_len = 30000
    every_nth = 2000

    train = data[:, shift:every_nth:shift+train_len-1]
    test = data[:, (shift+train_len):(shift+train_len+predict_len-1)]

    for rs in res_sizes
        for a in alphas
            for d in degrees
                for s in sigmas
                    for b in betas
                        for r in radii
                            create_and_test_esn(rs,r,d,tanh,s,b,a,NLADefault(),false, train, test, predict_len; make_plot=true)
                        end
                    end
                end
            end
        end
    end

end

#@time main()


res_size = 1000
radius = .7
degree = 100
||||||| merged common ancestors
rfp = "./test/input/reactions_postNN.csv"
sfp = "./test/input/species_postNN.csv"
icfp = "./test/input/initcond1.csv"

T=10. 
zeta = 1.
omega = 0.5
F_UV=1.
A_v=2.
E = 0.5
dens = 1e4
p = UCLCHEM.Parameters(zeta, omega, T, F_UV, A_v, E, dens)

#tspan = (0., 10^7 * 365. * 24. * 3600.)
tspan = (0., 10^2 * 365. * 24. * 3600.)

nw_prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan, rate_factor=1000)
prob = ODEProblem(nw_prob.network, nw_prob.u0, tspan)
sol = solve(prob, CVODE_BDF(), saveat = 24*360000*365)
sol2 = UCLCHEM.solve(nw_prob, dt=10000)
println("Problem solved with CVODE")

v = sol2.u
data = Matrix(hcat(v...))
shift = 1
train_len = 400000
predict_len = 1000
train = data[:, shift:shift+train_len-1]
test = data[:, shift+train_len:shift+train_len+predict_len-1]

approx_res_size = 500
radius = 2.8
degree = 200
=======
    T=10. 
    zeta = 1.
    omega = 0.5
    F_UV=1.
    A_v=2.
    E = 0.5
    dens = 1e4
    p = UCLCHEM.Parameters(zeta, omega, T, F_UV, A_v, E, dens)

    #tspan = (0., 10^7 * 365. * 24. * 3600.)
    tspan = (0., 10^6 * 365. * 24. * 3600.)

    nw_prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan)
    #prob = ODEProblem(nw_prob.network, nw_prob.u0, tspan)
    #sol = solve(prob, CVODE_BDF(), saveat = 24*360000*365)
    sol2 = UCLCHEM.solve(nw_prob, time_factor_pre_1000_years=20)
    println("Problem solved with CVODE")

    v = sol2.u
    data = Matrix(hcat(v...))
    shift = 1
    train_len = 100000
    predict_len = 30000
    every_nth = 2000

    train = data[:, shift:every_nth:shift+train_len-1]
    test = data[:, (shift+train_len):(shift+train_len+predict_len-1)]

    for rs in res_sizes
        for a in alphas
            for d in degrees
                for s in sigmas
                    for b in betas
                        for r in radii
                            create_and_test_esn(rs,r,d,tanh,s,b,a,NLADefault(),false, train, test, predict_len; make_plot=true)
                        end
                    end
                end
            end
        end
    end

end

#@time main()


res_size = 3013
radius = .7
degree = 100
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77
activation = tanh
<<<<<<< HEAD
alpha = 1.
sigma = .1
||||||| merged common ancestors
sigma = 0.7
beta = 0.0
alpha = 0.05 # leaking factor
=======
alpha = .70
sigma = 1.6
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77
nla_type = NLADefault()
extended_states = false
beta = 0.000001

#train = read_in_train_data("./test/data/solution1.csv")

#train = log10.(read_in_train_data("./test/data/solution2.csv", 100))
#test = log10.(read_in_train_data("./test/data/solution2.csv", 105))
#starting_point = log10.(read_first_rows("./test/data/solution2.csv", 100))

<<<<<<< HEAD
# USE THIS
train = read_in_train_data("./test/data/solution2.csv", 1000)
test = read_in_train_data("./test/data/solution2.csv", 1000)
starting_point = read_first_rows("./test/data/solution2.csv", 100)
||||||| merged common ancestors
@time W_out = ESNtrain(esn, beta)
@time output = ESNpredict(esn, predict_len, W_out)
=======
# USE THIS
#train = read_in_train_data("./test/data/solution2.csv", 1000)
#test = read_in_train_data("./test/data/solution2.csv", 1000)
#starting_point = read_first_rows("./test/data/solution2.csv", 100)
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77

#train = log10.(read_in_train_data("./test/data/solution2.csv", 1000))
#test = log10.(read_in_train_data("./test/data/solution2.csv", 1005))
#starting_point = log10.(read_first_rows("./test/data/solution2.csv", 100))

<<<<<<< HEAD
#train = read_in_train_data("./test/data/solution3.csv", 100)
#test = read_in_train_data("./test/data/solution3.csv", 105)
#starting_point = read_first_rows("./test/data/solution3.csv", 100)

train_subset = train[:, 2:2:998]
test_subset = train[:, 3:2:1000]

train_subset = train[:, 1:2000]
test_subset = train[:, 2001:end]
esn = ESN(res_size,
          train_subset,
          degree,
          radius,
          activation = activation, #default = tanh
          alpha = alpha, #default = 1.0
          sigma = sigma, #default = 0.1
          nla_type = nla_type, #default = NLADefault()
          extended_states = extended_states)


@time W_out = ESNtrain(esn, beta)
# reset the states
#esn.states[:, end] = esn.states[:, 2]
fake_state = zeros(Float64, esn.res_size, 1)
y = test_subset[:,1]
esn.states[:,end] =  (1-alpha).*fake_state + alpha*activation.((esn.W*fake_state)+(esn.W_in*y))
@time output = ESNpredict(esn, size(test_subset, 2), W_out)

plot(output[1:20,:]',layout=(4,5), label="predicted",size=(1600,900))
#plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
plot!(test_subset[1:20,2:end]',layout=(4,5), label="actual", size=(1600,900))
xaxis!(:log10)
||||||| merged common ancestors
using Plots
plot(transpose(output[1:12,:]),layout=(4,3), label="predicted",size=(1200,800))
plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
=======
#train = read_in_train_data("./test/data/solution3.csv", 100)
#test = read_in_train_data("./test/data/solution3.csv", 105)
#starting_point = read_first_rows("./test/data/solution3.csv", 100)

train_subset = train[:, 2:2:998]
test_subset = train[:, 3:2:1000]
esn = ESN(res_size,
          train_subset,
          degree,
          radius,
          activation = activation, #default = tanh
          alpha = alpha, #default = 1.0
          sigma = sigma, #default = 0.1
          nla_type = nla_type, #default = NLADefault()
          extended_states = extended_states)


@time W_out = ESNtrain(esn, beta)
# reset the states
esn.states[:, end] = esn.states[:, 2]
fake_state = zeros(Float64, esn.res_size, 1)
y = starting_point[:,3]
esn.states[:,end] =  (1-alpha).*fake_state + alpha*activation.((esn.W*fake_state)+(esn.W_in*y))
@time output = ESNpredict(esn, size(train, 2), W_out)

scatter(transpose(output[1:20,:]),layout=(4,5), label="predicted",size=(1600,900))
#plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
plot!(transpose(test_subset[1:20,2:end]),layout=(4,5), label="actual", size=(1600,900))
xaxis!(:log10)
>>>>>>> bb4e7ae053f485fcee256a5062df8875d01efe77
