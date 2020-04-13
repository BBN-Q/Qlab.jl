using Test, RandomQuantum, QuantumInfo, Cliffords

import QuantumInfo.eye

function check_LSStomo_result(true_values::Dict{String, Float64}, result::QuantumTomography.LSStateTomo)
    @test Set(keys(true_values)) == Set(keys(result.fit_params))
    @test Set(keys(true_values)) == Set(keys(result.errors))

    for k in keys(true_values)
        println("$(k): $(result.fit_params[k]) => $(true_values[k]) ± $(result.errors[k])")
        @test in_range(result.fit_params[k], true_values[k]-10*result.errors[k], true_values[k]+10*result.errors[k])
    end
end

function check_MLtomo_result(true_values::Dict{String, Float64}, result::QuantumTomography.MLStateTomo)
    @test Set(keys(true_values)) == Set(keys(result.fit_params))
    @test Set(keys(true_values)) == Set(keys(result.errors))

    for k in keys(true_values)
        println("$(k): $(result.fit_params[k]) => $(true_values[k]) ± $(result.errors[k])")
        @test in_range(result.fit_params[k], true_values[k]-10*result.errors[k], true_values[k]+10*result.errors[k])
    end
end

# setup single-qubit LSQ_QST and ML_QST
"""
Create a random single-qubit state along with its obserables
"""
function qst_1q_test_setup():
    # observe along the ⨦X, ⨦Z, ⨦Y basis
    obs = [ (complex(Pauli(i))+QuantumInfo.eye(2))/2 for i in 1:3 ]
    append!(obs, [ (-complex(Pauli(i))+QuantumInfo.eye(2))/2 for i in 1:3 ])

    ρ = .98*projector(rand(FubiniStudyPureState(2)))+.02*rand(HilbertSchmidtMixedState(2))

    return ρ, obs
end

"""
Create Auspex data from the random qubit state info
"""
function create_1q_data(ρ, obs):

end

function do_1q_QST( data):
    rhoLSQ,rhoML  = analyzeStateTomo(data,1,6)

end

function test_1q_QST():
    ρ, obs = qst_1q_test_setup()
    tomo = LSStateTomo(obs
    ideal_means = predict(tomo, ρ) |> real
    results = fit(tomo, ideal_means, ones(6))
end

function test_1q_ML_QST(; n=10_000, β=0.0, asymptotic=false, alt=false, maxiter=5000, ϵ=1000):
    ρ, obs = qst_1q_test_setup()
    tomo = MLStateTomo(obs,β)

    asymptotic_means = real(predict(tomo,ρ))

    samples = asymptotic ? asymptotic_means[1:3] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in asymptotic_means[1:3]]
    append!(samples, 1 .- samples)

    ρest, obj, status = fit(tomo, samples, maxiter=maxiter, δ=1/ϵ, λ=β)

    return status, trnorm(ρ-ρest), obj, ρest

end

# setup and test two-qubit LSQ_QST and ML_QST
"""
Create a random two-qubit state along with its obserables
"""
function qst_2q_test_setup()
    first_axes = Matrix[ (complex(kron(Pauli(i),Pauli(j)))+eye(4))/2 for j in 1:3, i in 1:3 ]
    obs = vcat(first_axes, Matrix[ (-complex(kron(Pauli(i), Pauli(j)))+eye(4))/2 for j in 1:3, i in 1:3 ])

    ρ = .98*projector(rand(FubiniStudyPureState(4)))+.02*rand(HilbertSchmidtMixedState(4))

    return ρ, obs
end

"""
Create Auspex data from the random two-qubit state info
"""
function create_2q_data(ρ, obs):

end

function test_2q_QST(data):
    rhoLSQ,rhoML  = analyzeStateTomo(data,2,6)
end

function test_2q_LSQ_acc(rhoLSQ, ρ, obs):
    ρ, obs = qst_2q_test_setup()
    tomo = LSStateTomo(obs)

end

function test_2q_ML_acc(rhoML, ρ, obs):
    ρ, obs = qst_2q_test_setup()
    ml_tomo = MLStateTomo(obs)

end
