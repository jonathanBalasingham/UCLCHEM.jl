#using Flux
using UCLCHEM
using Sundials
using Dates
using OrdinaryDiffEq
using ReservoirComputing
using ParameterizedFunctions
using Plots
using SQLite

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
                            predict_len)

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
    #p = plot(transpose(output[1:12,:]),layout=(4,3), label="predicted",size=(1200,800))
    #plot!(transpose(test[1:12,:]),layout=(4,3), label="actual", size=(1200,800))
    file_name = "ESN_$res_size" * "_$radius" * "_$degree" * "_$sigma" * "_$beta" * "_$alpha" * "_.png"
    print(file_name)
    #savefig(p, "./test/CTESN_Plots/" * file_name)
end


function main()
    res_sizes = [2999]
    radii = collect(.1:.1:10)
    degrees = [10, 50, 100]
    sigmas = collect(2.0:.2:20.0)
    betas = collect(0.0:.2:1.0)
    alphas = collect(.49:.1:.99)


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
    tspan = (0., 10^6 * 365. * 24. * 3600.)

    nw_prob = UCLCHEM.formulate(sfp,rfp,icfp,p,tspan, rate_factor=10000)
    #prob = ODEProblem(nw_prob.network, nw_prob.u0, tspan)
    #sol = solve(prob, CVODE_BDF(), saveat = 24*360000*365)
    sol2 = UCLCHEM.solve(nw_prob, time_factor_pre_1000_years=20)
    println("Problem solved with CVODE")

    v = sol2.u
    data = Matrix(hcat(v...))
    shift = 100000
    train_len = 3000000
    predict_len = 200000
    every_nth = 2000

    train = data[:, shift:every_nth:shift+train_len-1]
    test = data[:, (shift+train_len):(shift+train_len+predict_len-1)]

    for rs in res_sizes
        for a in alphas
            for d in degrees
                for s in sigmas
                    for b in betas
                        for r in radii
                            create_and_test_esn(rs,r,d,tanh,s,b,a,NLADefault(),false, train, test, predict_len)
                        end
                    end
                end
            end
        end
    end

end

@time main()
