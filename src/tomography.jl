using QuantumTomography, Cliffords, LinearAlgebra

"""
    tomo_gate_set(nbrQubits, nbrPulses; pulse_type, prep_meas)

Return a set of state preparation or readout unitary gates.
pulse_type: prepared states/meas. axes
prep_meas: 1 for prep., 2 for meas. pulses
"""
function tomo_gate_set(nbrQubits, nbrPulses; pulse_type="Clifford", prep_meas = 1)
    if nbrPulses==4
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
            error("Invalid prep./meas. pulse pulse_type")
        end
    elseif nbrPulses==6
        # Six pulse set
        Uset1Q = [complex(RI),
                  exp(-im*pi/4*X),
                  exp(+im*pi/4*X),
                  exp(-im*pi/4*Y),
                  exp(+im*pi/4*Y),
                  -im*X]
    elseif nbrPulses==12
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
        error("Invalid number of pulses");
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
    QST_LSQ(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps, n)

Function to perform least-squares inversion of state tomography data.

   + expResults : data array
   + varmat : covariance matrix for data
   + measPulseMap : array mapping each experiment to a measurement readout
     pulse
   + measOpMap: array mapping each experiment to a measurement channel
   + measPulseUs : array of unitaries of measurement pulses
   + measOps : array of measurement operators for each channel
"""
function QST_LSQ(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps)
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
    QST_ML(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps, n)

Function to perform maximum-likelihood quantum state tomography.

   + expResults : data array
   + varmat : covariance matrix for data
   + measPulseMap : array mapping each experiment to a measurement readout
     pulse
   + measOpMap: array mapping each experiment to a measurement channel
   + measPulseUs : array of unitaries of measurement pulses
   + measOps : array of measurement operators for each channel
"""
function QST_ML(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps)
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

function analyzeStateTomo(data::Dict{String,Dict{String,Array{Any,N} where N}}, nbrQubits, nbrPulses, nbrCalRepeats=2)

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
    measOpMap = repeat(1:numMeas, inner=nbrPulses^nbrQubits)
    measPulseMap = repeat(1:nbrPulses^nbrQubits, outer=numMeas)
    # Use a helper to get the measurement unitaries.
    measPulseUs = tomo_gate_set(nbrQubits, nbrPulses)

    # Now call the inversion routines
    # First least squares
    rhoLSQ = QST_LSQ(tomoData, varData, measPulseMap, measOpMap, measPulseUs, measOps)

    # Now constrained maximum-likelihood
    rhoML = QST_ML(tomoData, varData, measPulseMap, measOpMap, measPulseUs, measOps)
    # plotting to be implemented    pauliSetPlot(rho2pauli(rhoLSQ), newplot)

    return rhoLSQ, rhoML
end

"""
    rho2pauli(rho)

Convert a density matrix to a Pauli set vector.
"""
function rho2pauli(ρ)
    n = round(Int, log2(size(ρ,1)))
    paulis = sort(allpaulis(n), by=weight, dims=1)
    paulivec = [real(tr(ρ * p)) for p in paulis]
    return paulivec, paulis
end
