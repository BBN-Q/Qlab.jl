using QuantumTomography, Cliffords, LinearAlgebra, StatsBase

"""
    squeeze(A::AbstractArray)

Drops singleton dimentions from array structures.
Helper function to be used with caution!

https://stackoverflow.com/questions/52505760/dropping-singleton-dimensions-in-julia
"""
function squeeze(A::AbstractArray)
    singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)
    return dropdims(A, dims=singleton_dims)
end

"""
    tomo_gate_set(nbrQubits, nbrAxes; pulse_type::String="Clifford", prep_meas::Integer=1)

Return a set of state preparation or readout unitary gates for a give
`nbrQubits` along a given `nbrAxes`.
See http://arxiv.org/abs/quant-ph/0308098v1 for more information

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
- `varmat::Array{Number}`: covariance array for the expResults
- `measPulseMap::`: array mapping each experiment outcome to the corresponding
                    measPulseUs
- `measOpMap::`: array mapping each experiment outcome to a measurement operator
- `measPulseUs::`: array of unitaries applied before measurement pulses
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
        obs_ct = U' * op * U
        obs_ct = (obs_ct + obs_ct')/2 # force to be Hermitian
        push!(obs, obs_ct)
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
- `expResults::Array{Number}`: data array.  This is the list of expectation
                               values for each of the measurements obserables.
- `varmat::Array{Number}`: variance matrix for data in expResults
- `measPulseMap::`: array mapping each experiment to a measurement readout
                    pulse
- `measOpMap::`: array mapping each experiment to a measurement channel
- `measPulseUs::`: array of unitaries of measurement pulses.  These are
                   pulses applied before measurement in the experimental data,
                   mapping the measurement to the correct axis
- `measOps::`: array of measurement operators for each channel.  In the case of
               ML tomography, these need to be text book projectors to diagonal states.

# Returns
- `ρest::Array{Complex{Float64},d}`: d dimensional estimate of the density
                                     matrix obtained by constrained
                                     least-squares.

# Examples
```julia-repl
julia> QST_LSQ(expResults, measPulseMap, measOpMap, measPulseUs, measOps)
2×2 Array{Complex{Float64},2}:
 0.295239+2.94858e-14im  -0.21073+0.3957im
 -0.21073-0.3957im       0.704761-3.64935e-15im
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
        obs_ct = U' * op * U
        obs_ct = (obs_ct + obs_ct')/2 # force to be Hermitian
        push!(obs, obs_ct)
    end
    tomo = LSStateTomo(obs)

    ρest, obj, status = fit(tomo, expResults, varMat)
    if status != :Optimal
        println("LSStateTomo fit return status: $status")
    end
    return ρest
end

"""
    QST_ML(expResults, measPulseMap, measOpMap, measPulseUs, measOps)

Function to perform maximum-likelihood quantum state tomography.  This function
is usually wrapped by the analyzeStateTomo function.

# Arguments
- `expResults::Array{Number}`: data array.  This is a list of expectation
                               values for each of the measurements obserables.
- `measPulseMap::`: array mapping each experiment to a measurement
                    readout pulse
- `measOpMap::`: array mapping each experiment to a measurement channel
- `measPulseUs::`: array of unitaries of measurement pulses.  These are
                   pulses applied before measurement in the experimental data,
                   mapping the measurement to the correct axis
- `measOps::`: array of measurement operators for each channel.  In the case of
               ML tomography, these need to be text book projectors to
               diagonal states.

# Returns
- `ρest::Array{Complex{Float64},d}`: d dimensional estimate of the density
                                     matrix obtained by maximum likelihood.

# Examples
```julia-repl
julia> QST_ML(expResults, measPulseMap, measOpMap, measPulseUs, measOps)
2×2 Array{Complex{Float64},2}:
 0.295239+2.94858e-14im  -0.21073+0.3957im
 -0.21073-0.3957im       0.704761-3.64935e-15im
```
"""
function QST_ML(expResults,
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
    if length(obs) < 6
        @warn("State observations do not form a POVM!")
    end

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
- `nbrAxes`: number of measurements.  Either 4 or 6.  12 is possible but left
             for the user to do manually.
- `nbrCalRepeats`: number of repeated calibration points per calibration state

# Returns
- `rhoLSQ` : a 2^nbrQubits x 2^nbrQubits complex density matrix reconstructed
            with least-squares
- `rhoML` : (optional) a 2 x 2 complex density matrix reconstructed
            with least-squares.  Note this is only supported with single qubit
            data and six axes measurements.

# Examples
julia> datapath = "/path/to/data/folder"
julia> data, desc = load_data(datapath,24,"200315",load_var=true);
julia> rhoLSQ  = Qlab.analyzeStateTomo(data[1],2,4);
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
        if nbrCalRepeats == 0
            # In this case, assume data is already calibrated and insert
            # standard projectors
            append!(measOps, real(Qlab._create_ml_POVM(nbrQubits)))
            append!(tomoData, real(data_ql[1:end]))
            append!(varData, real(data_q["variance"])[1:end])
        else
            calData = real(data_ql[end-nbrCalRepeats*(2^nbrQubits )+1:end])
            avgCalData = mean(reshape(calData, nbrCalRepeats, 2^nbrQubits), dims=1)
            # Pull out the calibrations as diagonal measurement operators
            push!(measOps, diagm(0 => avgCalData[:]))

            #The data to invert
            append!(tomoData, real(data_ql[1:end-nbrCalRepeats*(2^nbrQubits)]))

            #variance
            append!(varData, real(data_q["variance"])[1:end-nbrCalRepeats*(2^nbrQubits)])
        end
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

    # Constrained maximum-likelihood is currently unsupported in general
    #
    # The reason for this has to do with the delicate nature of constructing
    # the correlations data and make operators like POVMs physical etc... If
    # this becomes a dire need, we can always revisit.  For most applications
    # of interest, LSQ is perfectly good

    # The one exception to this is single qubit data where measurements form a
    # proper POVM.  In experimental language this is the case where we take six
    # data points for the state reconstruction
    if nbrAxes == 6 && nbrQubits == 1
        rhoML = QST_ML(tomoData,
                       measPulseMap,
                       measOpMap,
                       measPulseUs,
                       measOps)

        return rhoLSQ, rhoML
    else
        return rhoLSQ, []
    end
    # plotting to be implemented    pauliSetPlot(rho2pauli(rhoLSQ), newplot)
end

"""
    _create_ml_POVM(numQubits::Int)

Create the textbook POVM for a given qubit system.

# Returns
- `mlPOVM::Array{ComplexF64}` : an array of numQubits^2 x numQubits^2
                                dimentional set of projectors that form a POVM

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
given tomography data set.  Helper function for the StateTomo structure.  This
function makes many assumptions about the structure of the data and is built
using heuristics.  YOur milage may vary

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

    # Assume zero cals as a base case for calibration or manual data scaling
    numCalRepeats = 0
    numCals = 0
    nbr_basis_states = (numQubits == 1) ? 2 : 4

    # Determine the cal repeats number
    # Given the number of qubits and the limited possible number of observables
    # search for a possible number of calibration repeats that matches the data
    for i in [4,6].^numQubits
        numCals_guess = numDataPoints - i
        nbrRepeats_guess = numCals_guess/nbr_basis_states
        if numDataPoints <= 6
            # catch the case where there are no cals
            nbrRepeats_guess = 0
        end
        if nbrRepeats_guess % 1 != 0
            # correct number will be a whole number
            continue
        end
        if nbrRepeats_guess < 0
            # filter negitive guesses
            continue
        end
        if !iseven(Int(nbrRepeats_guess))
            # This is a very safe assumption
            @warn("Assuming numCalRepeats is even!")
            continue
        end
        if !ispow2(Int(nbrRepeats_guess))
            # This is likely the case but warn the user
            @warn("Assuming nbr repeats is a power of 2!")
            continue
        end
        numCalRepeats = nbrRepeats_guess
    end
    if numQubits == 1 && numCalRepeats != 0
        numCals = numCalRepeats * 2
    elseif numQubits == 2 && numCalRepeats != 0
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

Determine the number of qubits, and return organized datasets based on the
structure of the data.  Single dimensional data is assumed to be averaged
and 2D data is assumed to be (experiment, shot) data.  This is designed to be
a private function of the StateTomo structure.  Note, any shot data found is
left to the user to process manually.
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
        ql = match(r"([qQ]\d?\d?\d?)[ -_](\w*)", i)
        push!(qubits, string(ql[1]))
        push!(labels, string(ql[2]))
    end
    # filter out any repeated entries
    unique!(qubits)
    unique!(labels)
    println("Found $(length(qubits)) sets of qubit data: " * string([string(i) * " " for i in qubits]...))
    println("Found $(length(labels)) sets of qubit data labels: " * string([string(i) * " " for i in labels]...))
    numQubits = length(qubits)
    numDatasets = length(labels)

    #pull out any correlation data
    correlation_data_sets = []
    qubit_correlation_keys = filter(x -> occursin(r"([Cc]orrelated?)", x), keys(data))
    if length(qubit_correlation_keys) != 0
        println("Correlation data found... ✓")
        corrData = true
    else
        println("Correlation data found... no")
        corrData = false
        if numQubits > 1
            @error("This appears to be two-qubit data but no correlation data
                    is provided!  Tomography will not work!  If you have
                    correlation data, please add it manually to the data
                    dictionary.")
        end
    end

    # check that at least one variance dataset exists
    # Least-squares reconstruction will not be possible without variance data
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

    # get the data size and classify assuming multi-dimensional data is raw,
    # integrated shot data
    tomo_data_idx = empty([], String)
    shot_data_idx = empty([], String)
    for i in keys(data)
        dims = size(squeeze(data[i]["data"]))
        if length(dims) == 1
            println("Main data set: " * string(i))
            push!(tomo_data_idx, i)
            numDataPoints = dims[1]
        elseif length(dims) > 1
            println("Shots data set: " * string(i))
            push!(shot_data_idx, i)
        end
    end

    tomoDataSets = filter(p -> p.first in tomo_data_idx, data)
    shotDataSets = filter(p -> p.first in shot_data_idx, data)

    return numQubits, numDatasets, corrData, varData, numDataPoints,
                                                      tomoDataSets,
                                                      shotDataSets
end
"""
State tomography object

This holds all the infromation necessary to do state and process tomography.

The object is constructed by passing it a tomography dataset and its descriptor
loaded from load_data.  Once created, this object can be passed directly to the
any of the tomographyic reconstruction methods.
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

"""
    analyzeStateTomo(tomo::StateTomo)

Function to setup and run quantum state tomography for a given number of qubits
and measurement axes.

# Arguments
- `tomo::Qlab.StateTomo`: State tomography object constructed from the data

# Returns
- `rhoLSQ`: a 2^nbrQubits x 2^nbrQubits complex density matrix reconstructed
            with least-squares
- `rhoML` : (optional) a 2 x 2 complex density matrix reconstructed
            with least-squares.  Note this is only supported with single qubit
            data and six axes measurements.

# Examples
julia> datapath = "/path/to/data/folder"
julia> data, desc = load_data(datapath,24,"200315",load_var=true);
julia> tomo = Qlab.StateTomo(data, desc);
julia> rhoLSQ  = Qlab.analyzeStateTomo(tomo);
julia> Qlab.pauli_set_plot(rhoLSQ)
"""
function analyzeStateTomo(tomoObj::StateTomo)

    data = tomoObj.tomoDataSets
    nbrQubits = tomoObj.numQubits
    nbrAxes = tomoObj.numAxes
    nbrCalRepeats = tomoObj.numCalRepeats
    measOps = Matrix{Float64}[]
    tomoData = Float64[]
    varData = Float64[]
    numMeas = length(data)

    # mlTomoData = Float64[]
    # mlVarData = Float64[]

    # For LSQ reconstruction - make the observables out of un-scaled
    # calibration points.  Note these do not form a POVM
    for data_q in values(data)
        # Average over calibration repeats
        data_ql = data_q["data"]
        if nbrCalRepeats == 0
            # In this case, assume data is already calibrated and insert
            # standard projectors
            append!(measOps, real(Qlab._create_ml_POVM(nbrQubits)))
            append!(tomoData, real(data_ql[1:end]))
            append!(varData, real(data_q["variance"][1:end]))
        else
            calData = real(data_ql[end-nbrCalRepeats*(2^nbrQubits )+1:end])
            avgCalData = mean(reshape(calData, nbrCalRepeats, 2^nbrQubits), dims=1)
            # Pull out the calibrations as diagonal measurement operators
            push!(measOps, diagm(0 => avgCalData[:]))

            #The data to invert
            append!(tomoData, real(data_ql[1:end-nbrCalRepeats*(2^nbrQubits)]))

            #variance
            append!(varData, real(data_q["variance"])[1:end-nbrCalRepeats*(2^nbrQubits)])
        end
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

    # Constrained maximum-likelihood is currently unsupported in general
    #
    # The reason for this has to do with the delicate nature of constructing
    # the correlation data and making operators, like POVMs, physical etc... If
    # this becomes a dire need, we can always revisit.  For most applications
    # of interest, LSQ is perfectly good.
    #
    # The one exception to this is single qubit data where measurements form a
    # proper POVM.  In experimental language this is the case where we take six
    # data points for the state reconstruction
    if nbrAxes == 6 && nbrQubits == 1
        rhoML = QST_ML(tomoData,
                       measPulseMap,
                       measOpMap,
                       measPulseUs,
                       measOps)

        return rhoLSQ, rhoML
    else
        return rhoLSQ, []
    end
end

"""
    rho2pauli(ρ)

Convert a density matrix, ρ, to a Pauli set vector of Pauli expectation values.

# Arguments
- `ρ`: State tomography object constructed from the data

# Returns
- `paulivec:Array{Float64}`: array of Pauli expectation vaules
- `paulis:Array{Pauli}`: array of Pauli operators as defined in Cliffords.jl

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
