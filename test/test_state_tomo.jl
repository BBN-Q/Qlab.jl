using Test,
      RandomQuantum,
      QuantumInfo,
      Cliffords,
      Distributions,
      QuantumTomography,
      Qlab

import SchattenNorms.trnorm

"""
Basic trace-norm comparison of the true state ρ and the estimate ρest
"""
function check_tomo_result(ρ::Array{<:Number,2},
                           ρest::Array{<:Number,2})
    @test isapprox(trnorm(ρ-ρest), 0, atol=3e-2)
end

################################################################################
# Setup tomography problem #####################################################
################################################################################

"""
Create the quantum state tomography problem for a given number of qubits

This setup defaults to the obervables expected in the python QGL package and
assumed in tht rest of the code base.

Returns a random density matrix, the set of observables, the mean values of
those observables and the variance in the case of least squares mean values.

If asymptotic is true, return the true predicted means.  If false, add noise.
"""
function qlab_qst_setup(num_Qubits;num_obs = 6,
                                   asymptotic=false,
                                   n = 10_000)

    if num_obs == 4
        # observe along the Z, X, -Y, -Z basis created in QGL
        obs_ = [(complex(Pauli(2)) + QuantumInfo.eye(2))/2,
               (complex(Pauli(3)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(1)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(2)) + QuantumInfo.eye(2))/2]
    elseif num_obs == 6
        # observe along the Z, X, -X, Y, -Y, -Z basis created in QGL
        obs_ = [(complex(Pauli(2)) + QuantumInfo.eye(2))/2,
               (complex(Pauli(3)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(3)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(1)) + QuantumInfo.eye(2))/2,
               (complex(Pauli(1)) + QuantumInfo.eye(2))/2,
               (complex(-Pauli(2)) + QuantumInfo.eye(2))/2]
    end

    # expand the state space based on the number of qubits
    obs = Array{ComplexF64,2}[]
    for gates in Base.product([obs_ for _ in 1:num_Qubits]...)
        push!(obs, kron(1, reverse(gates)...))
    end

    # create a random num_Qubit state
    ρ = .98*projector(rand(FubiniStudyPureState(2^num_Qubits))) +
                0.02*rand(HilbertSchmidtMixedState(2^num_Qubits))

    # restrict to least-squares if we only have four observables
    if num_obs == 4
        LSQtomo = LSStateTomo(obs)

        LSQ_means = real(predict(LSQtomo, ρ))

        # if asymptotic is true, add Binomial noise
        LSQ_samples = asymptotic ? LSQ_means[1:num_obs^num_Qubits] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in LSQ_means[1:num_obs^num_Qubits]]
        LSQ_var = asymptotic ? ones(length(LSQ_samples)) : n*(LSQ_samples - LSQ_samples.^2)/(n-1)

        return ρ, obs, LSQ_samples, LSQ_var, []

    elseif num_obs == 6
        LSQtomo = LSStateTomo(obs)
        MLtomo  = MLStateTomo(obs)

        LSQ_means = real(predict(LSQtomo, ρ))
        ML_means = real(predict(MLtomo, ρ))

        # if asymptotic is true, add Binomial noise
        LSQ_samples = asymptotic ? LSQ_means[1:num_obs^num_Qubits] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in LSQ_means[1:num_obs^num_Qubits]]
        LSQ_var = asymptotic ? ones(length(LSQ_samples)) : n*(LSQ_samples - LSQ_samples.^2)/(n-1)

        ML_samples = asymptotic ? ML_means[1:num_obs^num_Qubits] : Float64[rand(Distributions.Binomial(n,μ))/n for μ in ML_means[1:num_obs^num_Qubits]]

        return ρ, obs, LSQ_samples, LSQ_var, ML_samples
    end
end

################################################################################
# Test single-qubit LSQ_QST and ML_QST #########################################
################################################################################
"""
Test single qubit least squares tomography with the manual version of the
analyzeStateTomo function
"""
function do_1q_LSStateTomo(num_obs)

    ρ, obs, LSQ_data, LSQ_var, _ = qlab_qst_setup(1,
                                                  num_obs = num_obs,
                                                  asymptotic=false);

    # did reconstruction work without Qlab.jl??
    LSQ_tomo = LSStateTomo(obs);
    ρest, obj, status = fit(LSQ_tomo, LSQ_data, LSQ_var);
    @test status == :Optimal
    @test trnorm(ρ-ρest) < 3e-2

    # pack the data as the analyzeStateTomo function expects
    data = Dict{String,Dict{String,Array{Any,N} where N}}("q1-main"=>Dict("data"=>LSQ_data, "variance"=>LSQ_var))

    rhoLSQ, _  = Qlab.analyzeStateTomo(data,1,num_obs,nbrCalRepeats=0)
    check_tomo_result(ρ, rhoLSQ)
end

"""
Test single qubit, maximum likelihood by manually creating the set of
observables, data and projector using the QST_ML function
"""
function do_1q_MLStateTomo()
    # Note here we only test with num_obs = 6 because the num_obs = 4 case does
    # not form a POVM
    ρ, obs, _, _, ML_data = qlab_qst_setup(1,
                                           num_obs = 6,
                                           asymptotic=false);

    numQubits = 1
    numAxes = 6
    numCals = 0
    numMeas = 1
    proj = Qlab._create_ml_POVM(numQubits)

    # did reconstruction work without Qlab.jl??
    ML_tomo = MLStateTomo(obs);
    ρest, obj, status = fit(ML_tomo, ML_data);
    @test status == :Optimal
    @test trnorm(ρ-ρest) < 3e-2

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
    check_tomo_result(ρ, rhoML)
end

"""
Test the creation of the StateTomo object with single qubit least squares.
This test reads in data, builds the object and passes it to the analyzeStateTomo
function.
"""
function test_tomo_obj(num_obs)
    ρ, obs, LSQ_data, LSQ_var, _ = qlab_qst_setup(1,
                                                  num_obs = num_obs,
                                                  asymptotic=false);

    # did reconstruction work without Qlab.jl??
    LSQ_tomo = LSStateTomo(obs);
    ρest, obj, status = fit(LSQ_tomo, LSQ_data, LSQ_var);
    @test status == :Optimal
    @test trnorm(ρ-ρest) < 5e-2

    # pack the data as the analyzeStateTomo function expects
    data = Dict{String,Dict{String,Array{Any,N} where N}}("q1-main"=>Dict("data"=>LSQ_data, "variance"=>LSQ_var))
    desc = Dict{String,Any}("q1-main"=>Any[0])

    tomo = Qlab.StateTomo(data, desc)
    rhoLSQ, _  = Qlab.analyzeStateTomo(tomo)
    check_tomo_result(ρ, rhoLSQ)
end

################################################################################
# Test two-qubit LSQ_QST and ML_QST #########################################
################################################################################

# Note: only LSQ is supported`

function do_2q_LSStateTomo(num_obs)

    ρ, obs, LSQ_data, LSQ_var, _ = qlab_qst_setup(2,num_obs = num_obs,
                                                    asymptotic=false);

    # did reconstruction work without Qlab.jl??
    LSQ_tomo = LSStateTomo(obs);
    ρest, obj, status = fit(LSQ_tomo, LSQ_data, LSQ_var);
    @test status == :Optimal
    @test trnorm(ρ-ρest) < 3e-2

    # pack the data as the analyzeStateTomo function expects
    data = Dict{String,Dict{String,Array{Any,N} where N}}("q1-main"=>Dict("data"=>LSQ_data, "variance"=>LSQ_var))
    desc = Dict{String,Any}("q1-main"=>Any[0])

    rhoLSQ, _  = Qlab.analyzeStateTomo(data,2,num_obs,nbrCalRepeats=0)
    check_tomo_result(ρ, rhoLSQ)
end

# execute the tests
do_1q_LSStateTomo(4)
do_1q_LSStateTomo(6)
do_1q_MLStateTomo()
test_tomo_obj(4)
test_tomo_obj(6)
do_2q_LSStateTomo(6)
