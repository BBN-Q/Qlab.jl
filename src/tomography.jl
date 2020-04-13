using QuantumTomography, Cliffords, LinearAlgebra

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

# Examples
julia> datapath = "/path/to/data/folder"
julia> data = load_data(datapath,24,"200315",load_var=true); #26
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
    measOpMap = repeat(1:numMeas, inner=nbrAxes^nbrQubits)
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
