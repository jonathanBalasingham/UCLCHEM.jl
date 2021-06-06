using Sundials


k1 = 0:.1:1.0
k2 = 0:500:3e3
k3 = 0:100000:1e7
tspan = (0.0,1e5)

function rober(du,u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    du[1] = -k₁*y₁+k₃*y₂*y₃
    du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
    du[3] =  k₂*y₂^2
    nothing
end

prob = ODEProblem(rober,[0.5,0.3,0.2],tspan,[0.04,3e7,1e4])
sol = solve(prob, CVODE_BDF())

res_size = 500
radius = .2
degree = 100
activation = tanh
alpha = 1.
sigma = .1
nla_type = NLADefault()
extended_states = false
beta = 0.000001

number_of_rates = 3
saved_points = size(sol.t, 1)
save_points = sol.t

solution_filepath = "./data/roberson_data.csv"
weights_filepath = "./data/roberson_weights.csv"
current_group2 = [0]
current_group1 = zeros(Int32, saved_points) .+ 1

x0 = [0.5, 0.3, 0.2]

for p1 in k1
    for p2 in k2
        for p3 in k3
            local p = [p1, p2, p3]
            local prob = ODEProblem(rober, x0, tspan, p)
            local sol = solve(prob, CVODE_BDF(), saveat=save_points)
            local data = Matrix(hcat(sol.u...))

            if size(sol.t, 1) != saved_points
                println("Skipping parameter set: $p")
                continue
            end

            local esn = ESN(res_size,
                            data,
                            degree,
                            radius,
                            activation = activation,
                            alpha = alpha, 
                            sigma = sigma, 
                            nla_type = nla_type,
                            extended_states = extended_states)


            @time local W_out = ESNtrain(esn, beta)

            open(solution_filepath, "a") do io
                #writedlm(io, [current_group1 sol.t data'], ',')
            end
            weights = reshape(W_out', (1, size(W_out, 2)*size(W_out, 1)))
            open(weights_filepath, "a") do io
                writedlm(io, [current_group2 p' weights], ',')
            end

            current_group1 .+= 1
            current_group2 .+= 1
        end
    end
end