using QuantumTomography, Cliffords, LinearAlgebra, StatsBase

"""
    tomo_gate_set(nbrQubits, nbrAxes; pulse_type::String="Clifford", prep_meas::Integer=1)

Return a set of state preparation or readout unitary gates for a give
`nbrQubits` along a given `nbrAxes`.

# Arguments
- `nbrAxes::Integer`: number of single-qubit axes ∈ [4,6,12]
- `pulse_type::String="Clifford"`: prepared states/meas. axes
- `prep_meas::Integer=1`: 1 for prep gates, 2 for meas. gates

# Examples
```julia-repl
julia> tomo_gate_set(2, 4)
16-element Array{Array{Complex{Float64},2},1}:
[1.0 - 0.0im -0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im; -0.0 + 0.0im 1.0 - 0.0im
 0.0 + 0.0im -0.0 + 0.0im; -0.0 + 0.0im 0.0 + 0.0im 1.0 - 0.0im -0.0 + 0.0im;
 0.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im 1.0 - 0.0im]
[-0.7071067811865477 + 0.0im 0.0 + 0.7071067811865475im ...
```
"""
function tomo_gate_set(nbrQubits, nbrAxes; pulse_type="Clifford", prep_meas = 1)
    if nbrAxes==4
        # Four pulse set
        if pulse_type == "Clifford"
            Uset1Q = [complex(RI),
                      exp(-im*pi/4*X),
                      exp(-im*pi/4*Y),
                      -im*X]
        elseif pulse_type == "Tetra"
            if prep_meas == 1
                Uset1Q = [complex(RI),
                          exp(-im*acos(-1/3)*X),
                          exp(-im*2pi/3*Z)*exp(-im*acos(-1/3)*X),
                          exp(+im*2pi/3*Z)*exp(-im*acos(-1/3)*X)]
            else
                Uset1Q = [complex(RI),
                          exp(-im*acos(-1/3)*X),
                          exp(+im*acos(-1/3)*X)*exp(+im*2pi/3*Z),
                          exp(+im*acos(-1/3)*X)*exp(-im*2pi/3*Z)]
            end
        else
            error("Invalid prep./meas. pulse pulse_type.
            Must be ∈ {\"Tetra\", \"Clifford\"}")
        end
    elseif nbrAxes==6
        # Six pulse set
        Uset1Q = [complex(RI),
                  exp(-im*pi/4*X),
                  exp(+im*pi/4*X),
                  exp(-im*pi/4*Y),
                  exp(+im*pi/4*Y),
                  -im*X]
    elseif nbrAxes==12
        # 12 pulse set
        Uset1Q = [complex(RI),
                  -im*X,
                  -im*Y,
                  -im*Z,
                  exp(-im*pi/3*(+X+Y-Z)/sqrt(3)),  #X+Y-Z 120
                  exp(-im*pi/3*(+X-Y+Z)/sqrt(3)),  #X-Y+Z 120
                  exp(-im*pi/3*(-X+Y+Z)/sqrt(3)),  #-X+Y+Z 120
                  exp(-im*pi/3*(-X-Y-Z)/sqrt(3)),  #X+Y+Z -120 (equivalent to -X-Y-Z 120)
                  exp(-im*pi/3*(+X+Y+Z)/sqrt(3)),   #X+Y+Z 120
                  exp(-im*pi/3*(-X+Y-Z)/sqrt(3)),  #X-Y+Z -120 (equivalent to -X+Y-Z 120)
                  exp(-im*pi/3*(+X-Y-Z)/sqrt(3)),  #-X+Y+Z -120 (equivalent to X-Y-Z 120)
                  exp(-im*pi/3*(-X-Y+Z)/sqrt(3))]  #X+Y-Z -120 (equivalent to -X-Y+Z 120)
    else
        error("Invalid number of pulses.  Must be ∈ [4,6,12]");
    end

    # Now the gate set is the cartesian product of the 1Q gate set over the
    # number of qubits. Unfornately, Julia's default product is anti-
    # lexicographic (first index is fastest), so we need to reverse the
    # gate order before taking the kronecker product.
    gateSet = Array{ComplexF64,2}[]
    for gates in Base.product([Uset1Q for _ in 1:nbrQubits]...)
        push!(gateSet, kron(1, reverse(gates)...))
    end
    return gateSet
end

"""
    QST_LSQ(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps)

Function to perform unconstrained least-squares inversion of state
tomography data.

# Arguments
- `expResults::Array{Number}`: data array
- `varmat::Array{Number}`: covariance matrix for data
- `measPulseMap::`: array mapping each experiment to a measurement readout
     pulse
- `measOpMap::`: array mapping each experiment to a measurement channel
- `measPulseUs::`: array of unitaries of measurement pulses
- `measOps::`: array of measurement operators for each channel

# Examples
```julia-repl
julia> QST_LSQ(2, 4)

```
"""
function QST_FLSQ(expResults,
                  varMat,
                  measPulseMap,
                  measOpMap,
                  measPulseUs,
                  measOps)
    # construct the vector of observables for each experiment
    obs = Matrix{Complex{Float64}}[]
    for ct in 1:length(expResults)
        U = measPulseUs[measPulseMap[ct]]
        op = measOps[measOpMap[ct]]
        push!(obs, Hermitian(U' * op * U)) # force to be Hermitian
    end
    #return measOps
    # in order to constrain the trace to unity, add an identity observerable
    # and a corresponding value to expResults
    push!(obs, Diagonal((1.0)*fill(I.λ, size(measOps[1])))) # this can be replaced with I(size(measOps[1])) for Julia>=1.3
    expResults2 = [expResults; 1]
    # corresponding variance chosen arbitrarily (it should be very small)
    varMat2 = [varMat; minimum(varMat)]
    tomo = FreeLSStateTomo(obs)

    ρest, obj, status = fit(tomo, expResults2, varMat2)
    if status != :Optimal
        println("FreeLSStateTomo fit return status: $status")
    end
    return ρest
end
"""
    QST_LSQ(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps)

Function to perform constrained least-squares quantum state tomography.

# Arguments
- `expResults::Array{Number}`: data array
- `varmat::Array{Number}`: covariance matrix for data
- `measPulseMap::`: array mapping each experiment to a measurement readout
     pulse
- `measOpMap::`: array mapping each experiment to a measurement channel
- `measPulseUs::`: array of unitaries of measurement pulses
- `measOps::`: array of measurement operators for each channel

# Examples
```julia-repl
julia> QST_ML(2, 4)

```
"""
function QST_LSQ(expResults,
                 varMat,
                 measPulseMap,
                 measOpMap,
                 measPulseUs,
                 measOps)
    # construct the vector of observables for each experiment
    obs = Matrix{Complex{Float64}}[]
    for ct in 1:length(expResults)
        U = measPulseUs[measPulseMap[ct]]
        op = measOps[measOpMap[ct]]
        push!(obs, Hermitian(U' * op * U)) # force to be Hermitian
    end
    tomo = LSStateTomo(obs)

    ρest, obj, status = fit(tomo, expResults, varMat)
    if status != :Optimal
        println("LSStateTomo fit return status: $status")
    end
    return ρest
end

"""
    QST_ML(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps)

Function to perform maximum-likelihood quantum state tomography.

# Arguments
- `expResults::Array{Number}`: data array
- `varmat::Array{Number}`: covariance matrix for data
- `measPulseMap::`: array mapping each experiment to a measurement readout
     pulse
- `measOpMap::`: array mapping each experiment to a measurement channel
- `measPulseUs::`: array of unitaries of measurement pulses
- `measOps::`: array of measurement operators for each channel

# Examples
```julia-repl
julia> QST_ML(2, 4)

```
"""
function QST_ML(expResults,
                varMat,
                measPulseMap,
                measOpMap,
                measPulseUs,
                measOps; n=100_000, β=0.0, maxiter=5000, ϵ=1000)
    # construct the vector of observables for each experiment
    obs = Matrix{Complex{Float64}}[]
    for ct in 1:length(expResults)
        U = measPulseUs[measPulseMap[ct]]
        op = measOps[measOpMap[ct]]
        push!(obs, Hermitian(U' * op * U)) # force to be Hermitian
    end
    # ML obs must be POVMs -> Hermitian, possitive-semidefinite, and trace 1
    tomo = MLStateTomo(obs)

    ρest, obj, status = fit(tomo,
                            expResults,
                            maxiter=maxiter,
                            δ=1/ϵ,
                            λ=β)
    if status != :Optimal
        println("MLStateTomo fit return status: $status")
    end
    return ρest
end

"""
    analyzeStateTomo(data::Dict{String,Dict{String,Array{Any,N} where N}},
                          nbrQubits::Int,
                          nbrAxes::Int ∈ [4,6,12];
                          nbrCalRepeats::Int=2)

Function to setup and run quantum state tomography for a given number of qubits
and measurement axes.

# Arguments
- `data`: data array strucured as a set of data with a string name with a 'data'
          key and a 'variance' key.  The variance is required for all tomography
          reconstructions except for the free-LSQ tomo.  Also, for two-qubit
          tomography, the correlation data between the two qubit data sets is
          required for reconstruction.
- `nbrQubits`: number of qubits
- `nbrAxes`: number of measurements
- `nbrCalRepeats`: number of repeated calibration points per calibration state

# Returns
- `rhoLSQ`: a 2^nbrQubits x 2^nbrQubits complex density matrix reconstructed
            with least-squares
- `rhoML`: a 2^nbrQubits x 2^nbrQubits complex density matrix reconstructed
            with maximun-likelyhood

# Examples
julia> datapath = "/path/to/data/folder"
julia> data, desc = load_data(datapath,24,"200315",load_var=true);
julia> rhoLSQ,rhoML  = Qlab.analyzeStateTomo(data[1],2,4);
julia> Qlab.pauli_set_plot(rhoLSQ)
"""
function analyzeStateTomo(data::Dict{String,Dict{String,Array{Any,N} where N}},
                          nbrQubits,
                          nbrAxes;
                          nbrCalRepeats=2)

    measOps = Matrix{Float64}[]
    tomoData = Float64[]
    varData = Float64[]
    numMeas = length(data)

    for data_q in values(data)
        # Average over calibration repeats
        data_ql = data_q["data"]
        calData = real(data_ql[end-nbrCalRepeats*(2^nbrQubits )+1:end])
        avgCalData = mean(reshape(calData, nbrCalRepeats, 2^nbrQubits), dims=1)
        # Pull out the calibrations as diagonal measurement operators
        push!(measOps, diagm(0 => avgCalData[:]))

        #The data to invert
        append!(tomoData, real(data_ql[1:end-nbrCalRepeats*(2^nbrQubits)]) )

        #variance
        append!(varData, real(data_q["variance"])[1:end-nbrCalRepeats*(2^nbrQubits)] )
    end
    # Map each experiment to the appropriate readout pulse
    # These are integers 1:length(data), each repeated numAxes^nbrQubits times
    measOpMap = repeat(1:numMeas, inner=nbrAxes^nbrQubits)
    # These are integers 1:nbrAxes^nbrQubits, unrolled length(data) times
    measPulseMap = repeat(1:nbrAxes^nbrQubits, outer=numMeas)
    # Use a helper to get the measurement unitaries.
    measPulseUs = tomo_gate_set(nbrQubits, nbrAxes)

    # Now call the inversion routines
    # First least squares
    rhoLSQ = QST_LSQ(tomoData,
                     varData,
                     measPulseMap,
                     measOpMap,
                     measPulseUs,
                     measOps)

    # Now constrained maximum-likelihood
    rhoML = QST_ML(tomoData,
                   varData,
                   measPulseMap,
                   measOpMap,
                   measPulseUs,
                   measOps)
    # plotting to be implemented    pauliSetPlot(rho2pauli(rhoLSQ), newplot)

    return rhoLSQ, rhoML
end

"""
    _create_ml_POVM(numQubits::Int)

Create the textbook POVM for a given qubit system
"""
function _create_ml_POVM(numQubits::Int)
    mlPOVM = Array{ComplexF64}[]
    if numQubits == 1
        push!(mlPOVM, diagm([1.,0.]))
        push!(mlPOVM, diagm([0.,1.]))
    elseif numQubits == 2
        push!(mlPOVM, diagm([1.,0.,0.,0.]))
        push!(mlPOVM, diagm([0.,1.,0.,0.]))
        push!(mlPOVM, diagm([0.,0.,1.,0.]))
        push!(mlPOVM, diagm([0.,0.,0.,1.]))
    end
    return mlPOVM
end

"""
    _parese_exp_num(n::Int, n_qubits::Int)

Determine the number of calibration points, calibration repeats, axes in a
given tomography data set.

# Arguments
- `n` : total number of experimental data points including
        calibration points
- `n_qubits` : number of qubits represented in the data set.  Only one
               and two qubit tomography is supported.
# Returns
- `numCalRepeats::Int` : number of repeats for each calibration point
- `numCals::Int` : the total number of calibration points in an experiment
- `numAxes::Int` : the number of of axes were observations were made.
                    This must be 4 or 6.
"""
function _parese_exp_num(numDataPoints::Int, numQubits::Int)
    numCalRepeats = 0
    nbr_basis_states = (numQubits == 1) ? 2 : 4
    # determine the cal repeats number
    for i in [4,6].^numQubits
        numCals_guess = numDataPoints - i
        nbrRepeats_guess = numCals_guess/nbr_basis_states
        if nbrRepeats_guess % 1 != 0
            # correct number will be a whole number
            continue
        end
        if !isposdef(nbrRepeats_guess)
            # filter out obviously wrong guesses
            continue
        end
        if !iseven(Int(nbrRepeats_guess))
            @warn("Assuming numCalRepeats is even!")
            continue
        end
        if !ispow2(Int(nbrRepeats_guess))
            @warn("Assuming nbr repeats is a power of 2!")
            continue
        end
        numCalRepeats = nbrRepeats_guess
    end
    if numQubits == 1
        numCals = numCalRepeats * 2
    elseif numQubits == 2
        numCals = numCalRepeats * 4
    end

    #determine the number of axes
    if numQubits == 2
        numAxes = sqrt(numDataPoints-numCals)
    elseif numQubits ==1
        numAxes = numDataPoints-numCals
    end
    # assert numAxes must equal [4,6]
    if !(numAxes in [4,6])
        error("Obervables must be 4 or 6.  Please check your data!")
    end
    return numCals, numCalRepeats, numAxes
end

"""
    _pre_process_data(data::Dict{String,Dict{String,Array{Any,N} where N}},
                           desc::Dict{String,Any})
Preprocess the data

determine the number of qubits, and return organized datasets based on the
structure of the data.  Single dimensional data is assumed to be averaged
and 2D data is assumed to be (experiment, shot) data
"""
function _pre_process_data(data::Dict{String,Dict{String,Array{Any,N} where N}},
                           desc::Dict{String,Any})
    numDataPoints = 0
    corrData = false
    varData = false

    # load  and parse the data
    println("Preprocessing data")
    qubit_data_keys = filter(x -> occursin(r"([qQ]\d?\d?\d?)[ -_](\w*)", x), keys(data))

    qubits = []
    labels = []
    for i in qubit_data_keys
        foo = match(r"([qQ]\d?\d?\d?)[ -_](\w*)", i)
        push!(qubits, string(foo[1]))
        push!(labels, string(foo[2]))
    end
    unique!(qubits)
    unique!(labels)
    println("Found $(length(qubits)) sets of qubit data: " * string([string(i) * " " for i in qubits]...))
    println("Found $(length(labels)) sets of qubit data labels: " * string([string(i) * " " for i in labels]...))
    numQubits = length(qubits)
    numDatasets = length(labels)

    #pull out any correlation
    correlation_data_sets = []
    qubit_correlation_keys = filter(x -> occursin(r"([Cc]orrelated?)", x), keys(data))
    if length(qubit_correlation_keys) != 0
        println("Correlation data found... ✓")
        corrData = true
    else
        println("Correlation data found... no")
        corrData = false
    end

    #check that atleast one variance dataset exists
    variance_data = []
    for i in keys(data)
        if length(filter(x -> occursin(r"([Vv]ariance)", x), keys(data[i]))) != 0
            println("Variance data found for dataset: $(i)")
            varData = true
        else
            println("Variance data for $(i) found... no")
            varData = false
            @warn("You may not be able to do state tomography")
        end
    end

    # get the data size and classify them
    tomo_data_idx = empty([], String)
    shot_data_idx = empty([], String)
    for i in keys(data)
        data_size = size(data[i]["data"])
        if data_size[1] == 1
            println("Main data set: " * string(i))
            push!(tomo_data_idx, i)
            numDataPoints = data_size[2]
        elseif data_size[1] > 1
            println("Shots data set: " * string(i))
            push!(shot_data_idx, i)
        end
    end
#     tomoDataSets = filter((k,v) -> k in tomo_data_idx, data)
#     shotDataSets = filter((k,v) -> k in shot_data_idx, data)
    tomoDataSets = filter(p -> p.first in tomo_data_idx, data)
    shotDataSets = filter(p -> p.first in shot_data_idx, data)
    return numQubits, numDatasets, corrData, varData, numDataPoints,
                                                      tomoDataSets,
                                                      shotDataSets
end
"""
State tomography object

This holds all the infromation necessary to do state and process tomography.

The object is constructed by passing it a tomography data file or a dataset
and its descriptor loaded from load_data.  Once created, this object can be
passed directly to the any of the tomographyic reconstruction methods.
"""
struct StateTomo
    numQubits::Int
    numDatasets::Int
    corrData::Bool
    varData::Bool

    tomoDataSets::Dict{String,Dict{String,Array{Any,N} where N}}
    shotDataSets::Dict{String,Dict{String,Array{Any,N} where N}}

    numAxes::Int
    numCals::Int
    numCalRepeats::Int
    numDataPoints::Int

    mlPOVM::Array{Matrix{ComplexF64}}

    """
    Basic constructor
    """
    function StateTomo(data::Dict{String,Dict{String,Array{Any,N} where N}},
                        desc::Dict{String,Any})

        numQubits,
        numDatasets,
        corrData,
        varData,
        numDataPoints,
        tomoDataSets,
        shotDataSets = _pre_process_data(data, desc)

        ############################################################
        numCals, numCalRepeats, numAxes = _parese_exp_num(numDataPoints,
                                                          numQubits)
        ################################################################
        mlPOVM = _create_ml_POVM(numQubits)
        ################################################################
        new(numQubits, numDatasets, corrData, varData, tomoDataSets,
                                                       shotDataSets,
                                                       numAxes,
                                                       numCals,
                                                       numCalRepeats,
                                                       numDataPoints,
                                                       mlPOVM)
    end
end

function analyzeStateTomo(tomoObj::StateTomo)

    data = tomoObj.tomoDataSets
    nbrQubits = tomoObj.numQubits
    nbrAxes = tomoObj.numAxes
    nbrCalRepeats = tomoObj.numCalRepeats
    measOps = Matrix{Float64}[]
    tomoData = Float64[]
    mlTomoData = Float64[]
    varData = Float64[]
    mlVarData = Float64[]
    numMeas = length(data)

    # For LSQ reconstruction - make the observables out of un-scaled
    # calibration points.  Note these do not form a POVM
    for data_q in values(data)
        # Average over calibration repeats
        data_ql = data_q["data"]
        calData = real(data_ql[end-nbrCalRepeats*(2^nbrQubits )+1:end])
        avgCalData = mean(reshape(calData, nbrCalRepeats, 2^nbrQubits), dims=1)
        # Pull out the calibrations as diagonal measurement operators
        push!(measOps, diagm(0 => avgCalData[:]))

        #The data to invert
        append!(tomoData, real(data_ql[1:end-nbrCalRepeats*(2^nbrQubits)]) )

        #variance
        append!(varData, real(data_q["variance"])[1:end-nbrCalRepeats*(2^nbrQubits)] )
    end
    # Map each experiment to the appropriate readout pulse
    # These are integers 1:length(data), each repeated numAxes^nbrQubits times
    measOpMap = repeat(1:numMeas, inner=nbrAxes^nbrQubits)
    # These are integers 1:nbrAxes^nbrQubits, unrolled length(data) times
    measPulseMap = repeat(1:nbrAxes^nbrQubits, outer=numMeas)
    # Use a helper to get the measurement unitaries.
    measPulseUs = Qlab.tomo_gate_set(nbrQubits, nbrAxes)

    # Now call the inversion routines
    # First least squares
    rhoLSQ = QST_LSQ(tomoData,
                     varData,
                     measPulseMap,
                     measOpMap,
                     measPulseUs,
                     measOps)

    # Now constrained maximum-likelihood
    # the data needs to be scaled between -1 and 1
    # the maixmum likelyhood POVM also needs to be used

    for (i, data_q) in enumerate(values(data))
        # Average over calibration repeats

        # check for multi qubit data
        # we'll have to assume a qubit ordering for now
        if i == 1 && nbrQubits > 2
            orderBit = 2
        else
            orderBit = 1
        end
        data_ql = Qlab.cal_data(data_q["data"], bit=orderBit,
                                                    nqubits=nbrQubits,
                                                    num_repeats=nbrCalRepeats)
        #println(length(data_ql))
        #The data to invert
        append!(mlTomoData, real(data_ql) )

    end

    rhoML = QST_ML(mlTomoData,
                   varData,
                   measPulseMap,
                   measOpMap,
                   measPulseUs,
                   tomoObj.mlPOVM)

    return rhoLSQ, rhoML
end

"""
    rho2pauli(ρ)

Convert a density matrix, ρ, to a Pauli set vector.

# Examples
```julia-repl
julia> using Cliffords,
julia> foo = randn(ComplexF64, 2, 2)
2×2 Array{Complex{Float64},2}:
  0.310878-1.27868im    0.28776+1.87167im
 -0.904544+0.015083im  0.272544+0.620223im
julia> rho2pauli(foo)
([0.5834214596710264, -0.616783378678339, -1.8565830050748344,
 0.03833362063824197], Pauli{1}[+I, +X, +Y, +Z])
```
"""
function rho2pauli(ρ)
    n = round(Int, log2(size(ρ,1)))
    if n == 1 # special case for vectors
        paulis = sort(allpaulis(n), by=weight)
    else
        paulis = sort(allpaulis(n), dims=1, by=weight)
    end
    paulivec = [real(tr(ρ * p)) for p in paulis]
    return paulivec, paulis
end
