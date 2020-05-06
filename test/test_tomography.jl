using Test, RandomQuantum, QuantumInfo, Cliffords, Distributions

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
function qst_1q_test_setup(;asymptotic=false, n = 10_000)
    # observe along the ⨦X, ⨦Z, ⨦Y basis
    obs = [ (complex(Pauli(i))+QuantumInfo.eye(2))/2 for i in 1:3 ]
    append!(obs, [ (-complex(Pauli(i))+QuantumInfo.eye(2))/2 for i in 1:3 ])

    ρ = .98*projector(rand(FubiniStudyPureState(2)))+0.02*rand(HilbertSchmidtMixedState(2))

    LSQtomo = LSStateTomo(obs)
    MLtomo  = MLStateTomo(obs)

    LSQ_means = real(predict(LSQtomo, ρ))
    ML_means = real(predict(MLtomo, ρ))

    LSQ_samples = asymptotic ? LSQ_means[1:3] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in LSQ_means[1:3]]
    append!(LSQ_samples, 1 .- LSQ_samples)
    LSQ_var = asymptotic ? ones(length(samples)) : n*(LSQ_samples - LSQ_samples.^2)/(n-1)

    ML_samples = asymptotic ? ML_means[1:3] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in ML_means[1:3]]
    append!(ML_samples, 1 .- ML_samples)

    return ρ, obs, LSQ_samples, LSQ_var, ML_samples
end

# now if we use the Qlab.jl code??
function qlab_qst_1q_test_setup(;num_obs = 6, asymptotic=false, n = 10_000)
    if num_obs == 4
        # observe along the Z, X, Y, -Z basis created in QGL
        obs = [(complex(Pauli(2)) + QuantumInfo.eye(2))/2,
               (complex(Pauli(3)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(1)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(2)) + QuantumInfo.eye(2))/2]
    elseif num_obs == 6
        # observe along the Z, X, -X, Y, -Y, -Z basis created in QGL
        obs = [(complex(Pauli(2)) + QuantumInfo.eye(2))/2,
               (complex(Pauli(3)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(3)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(1)) + QuantumInfo.eye(2))/2,
               (complex(Pauli(1)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(2)) + QuantumInfo.eye(2))/2]
    end

    ρ = .98*projector(rand(FubiniStudyPureState(2))) +
                0.02*rand(HilbertSchmidtMixedState(2))

    LSQtomo = LSStateTomo(obs)
    MLtomo  = MLStateTomo(obs)

    LSQ_means = real(predict(LSQtomo, ρ))
    ML_means = real(predict(MLtomo, ρ))

    LSQ_samples = asymptotic ? LSQ_means[1:num_obs] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in LSQ_means[1:num_obs]]
    LSQ_var = asymptotic ? ones(length(LSQ_samples)) : n*(LSQ_samples - LSQ_samples.^2)/(n-1)

    ML_samples = asymptotic ? ML_means[1:num_obs] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in ML_means[1:num_obs]]

    return ρ, obs, LSQ_samples, LSQ_var, ML_samples
end

function do_1q_LSStateTomo(data, num_obs):
    # did reconstruction work without Qlab.jl??
    LSQ_tomo = LSStateTomo(obs);
    ρest, obj, status = fit(LSQ_tomo, LSQ_data, LSQ_var);

    ρ, obs, LSQ_data, LSQ_var, _ = qlab_qst_1q_test_setup(num_obs = num_obs,
                                                          asymptotic=false);
    # pack the data as the analyzeStateTomo function expects
    data = Dict("q1-main"=>Dict("data"=>LSQ_data, "variance"=>LSQ_var))
    rhoLSQ, _  = analyzeStateTomo(data,1,num_obs)
    return trnorm(ρ-rhoLSQ)
end

function do_1q_MLStateTomo(data, num_obs):
    # did reconstruction work without Qlab.jl??
    ML_tomo = MLStateTomo(obs);
    ρest, obj, status = fit(ML_tomo, ML_data);
    
    ρ, obs, _, _, ML_data = qlab_qst_1q_test_setup(num_obs = num_obs,
                                                          asymptotic=false);
    numQubits = 1
    numAxes = num_obs
    numCals = 0
    numMeas = 1
    proj = Qlab._create_ml_POVM(numQubits)

    # Map each experiment to the appropriate readout pulse
    # These are integers 1:length(data), each repeated numAxes^nbrQubits times
    measOpMap = repeat(1:numMeas, inner=numAxes^numQubits)
    # These are integers 1:nbrAxes^nbrQubits, unrolled length(data) times
    measPulseMap = repeat(1:numAxes^numQubits, outer=numMeas)
    # Use a helper to get the measurement unitaries.
    measPulseUs = Qlab.tomo_gate_set(numQubits, numAxes);

    rhoML = Qlab.QST_ML(ML_data,
                     measPulseMap,
                     measOpMap,
                     measPulseUs,
                     proj)
    return trnorm(ρ-rhoML)
end

# setup and test two-qubit LSQ_QST
# Note: only LSQ is supported`
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
function create_2q_data(ρ, obs)

end

function test_2q_QST(data):
    rhoLSQ,rhoML  = analyzeStateTomo(data,2,6)
end

function test_2q_LSQ_acc(rhoLSQ, ρ, obs)
    ρ, obs = qst_2q_test_setup()
    tomo = LSStateTomo(obs)

end

function test_2q_ML_acc(rhoML, ρ, obs)
    ρ, obs = qst_2q_test_setup()
    ml_tomo = MLStateTomo(obs)

end
