using SCS, QuantumInfo, Convex, Cliffords

"""
    tomo_gate_set(nbrQubits, nbrPulses; pulse_type, prep_meas)

Return a set of state preparation or readout unitary gates.
pulse_type: prepared states/meas. axes
prep_meas: 1 for prep., 2 for meas. pulses
"""
function tomo_gate_set(nbrQubits, nbrPulses; pulse_type="Clifford", prep_meas = 1)
  id, sx, sy, sz = Matrix[complex(allpaulis(1)[i]) for i in 1:4]
  gateSet = Dict()
  Uset1Q = Dict{Int, Matrix{Complex128}}()
  if nbrPulses==4
    #Four pulse set
    if pulse_type == "Clifford"
      Uset1Q[1] = id
      Uset1Q[2] = expm(-1im*(pi/4)*sx)
      Uset1Q[3] = expm(-1im*(pi/4)*sy)
      Uset1Q[4] = expm(-1im*(pi/2)*sx)
    elseif pulse_type == "Tetra"
      if prep_meas==1
        Uset1Q[1] = id
        Uset1Q[2] = expm(-1im*acos(-1/3)*sx)
        Uset1Q[3] = expm(-1im*2*pi/3*sz)*expm(-1im*acos(-1/3)*sx)
        Uset1Q[4] = expm(1im*2*pi/3*sz)*expm(-1im*acos(-1/3)*sx)
      else
        Uset1Q[1] = id
        Uset1Q[2] = expm(1im*acos(-1/3)*sx)
        Uset1Q[3] = expm(1im*acos(-1/3)*sx)*expm(1im*2*pi/3*sz)
        Uset1Q[4] = expm(1im*acos(-1/3)*sx)*expm(-1im*2*pi/3*sz)
      end
    else
      error("Invalid prep./meas. pulse pulse_type")
    end
  elseif nbrPulses==6
    #Six pulse set
    Uset1Q[1] = id
    Uset1Q[2] = expm(-1im*(pi/4)*sx)
    Uset1Q[3] = expm(1im*(pi/4)*sx)
    Uset1Q[4] = expm(-1im*(pi/4)*sy)
    Uset1Q[5] = expm(1im*(pi/4)*sy)
    Uset1Q[6] = expm(-1im*(pi/2)*sx)
  elseif nbrPulses==12
    #12 pulse set
    Uset1Q[1] = id
    Uset1Q[1] = expm(-1im*(pi/2)*sx)
    Uset1Q[1] = expm(-1im*(pi/2)*sy)
    Uset1Q[1] = expm(-1im*(pi/2)*sz)
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(sx+sy-sz))  #sx+sy-sz 120
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(sx-sy+sz))  #sx-sy+sz 120
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(-sx+sy+sz))  #-sx+sy+sz 120
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(-sx-sy-sz))  #sx+sy+sz -120 (equivalent to -sx-sy-sz 120)
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(sx+sy+sz))   #sx+sy+sz 120
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(-sx+sy-sz))  #sx-sy+sz -120 (equivalent to -sx+sy-sz 120)
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(sx-sy-sz))  #-sx+sy+sz -120 (equivalent to sx-sy-sz 120
    Uset1Q[1] = expm(-1im*(pi/3)*(1/sqrt(3))*(-sx-sy+sz))  #sx+sy-sz -120 (equivalent to -sx-sy+sz 120
  else
    error("Invalid number of pulses");
end

#Now kron things together
#First create a matrix with giving the mod nbrPulses description of which 1Q gates to kron together
numGates = nbrPulses^nbrQubits;
kronMat = zeros(UInt8, numGates, nbrQubits);
for qubitct = 1:nbrQubits
  kronMat[:,qubitct] = reshape(repmat(transpose(1:nbrPulses), nbrPulses^(nbrQubits-qubitct), nbrPulses^(qubitct-1)), numGates, 1);
end

for gatect = 1:numGates
  gateSet[gatect] = 1;
  for qubitct = 1:nbrQubits
    gateSet[gatect] = kron(gateSet[gatect], Uset1Q[kronMat[gatect, qubitct]]);
  end
end
    return gateSet
end

function QST_LSQ(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps, n)

#Function to perform least-squares inversion of state tomography data
#
# expResults : data array
# varmat : convariance matrix for data
# measPulseMap: array mapping each experiment to a measurement readout
# pulse
# measOpMap: array mapping each experiment to a measurement channel
# measPulseUs : cell array of unitaries of measurement pulses
# measOps : cell array of measurement operators for each channel
# n : number of qubits

#Construct the predictor matrix.  Each row is an experiment.  The number of
#columns is 4^n for the size of the vectorized density matrix.
#%First transform the measurement operators by the readout pulses to create
#the effective measurement operators and then flatten into row of the
#predictor matrix
predictorMat = complex(zeros(length(expResults), 4^n));
for expct = 1:length(expResults)

        tmp = transpose(measPulseUs[Int(measPulseMap[Int(expct)])]'*measOps[Int(measOpMap[Int(expct)])]*measPulseUs[Int(measPulseMap[Int(expct)])]);
        predictorMat[expct,:] = tmp[:];
end
    invVarMat = inv(varMat);
    A = predictorMat'*invVarMat*predictorMat
    B = predictorMat'*invVarMat*expResults
    rhoLSQ = A\B
    rhoLSQ = reshape(rhoLSQ, 2^n, 2^n);
    return rhoLSQ
end

function QST_SDP(expResults, varMat, measPulseMap, measOpMap, measPulseUs, measOps, n)
#Function to perform constrained SDP optimization of a physical density matrix
#consitent with the data.
#
# expResults : structure array (length total number of experiments)
#   each structure containts fields data, measPulse, measOperator
# measPulses : cell array of unitaries of measurment pulses
# measOps : cell array of measurment operators
# n : number of qubits

#Construct the predictor matrix.  Each row is an experiment.  The number of
#columns is 4^n for the size of the vectorized density matrix.
#First transform the measurement operators by the readout pulses to create
#the effective measurement operators and then flatten into row of the
#predictor matrix

 predictorMat = complex(zeros(length(expResults), 4^n));for expct = 1:length(expResults)
  tmp = transpose(measPulseUs[Int(measPulseMap[Int(expct)])]'*measOps[Int(measOpMap[Int(expct)])]*measPulseUs[Int(measPulseMap[Int(expct)])]);
     predictorMat[expct,:] = tmp[:];
end

solver = SCSSolver(verbose=0, max_iters=10_000, eps = 1e-8)

invVarMat = inv(varMat);
invVarMat = real(sqrtm(real(sqrtm(invVarMat'*invVarMat))));


ρr = Variable(2^n, 2^n)
ρi = Variable(2^n, 2^n)

constraints = trace(ρr) == 1
constraints += trace(ρi) == 0
constraints += isposdef([ρr ρi; -ρi ρr])

#We want to minimize the difference between predicted results and experimental results
problem = minimize( vecnorm(invVarMat*(expResults - [real(predictorMat) imag(predictorMat)]*[vec(ρr); vec(ρi)]), 2)^2, constraints )
solve!(problem, solver)

@printf("Done\n")
return (ρr.value - 1im*ρi.value)
end


function analyzeStateTomo(data, nbrQubits, nbrPulses, nbrCalRepeats)

  numMeas = length(data)
  measOps = Dict()
  tomoData = Array(Float64[])
  varData = Array(Float64[])
  if isa(data, Dict)
    datatemp = Dict()
    datatemp[1] = data
    numMeas = 1
  else
    datatemp = data
  end

  for ct = 1:numMeas
    #Average over calibration repeats
    calData = mean(reshape(real(datatemp[ct]["data"][end-nbrCalRepeats*(2^nbrQubits)+1:end]), nbrCalRepeats, 2^nbrQubits), 1);
    #Pull out the calibrations as diagonal measurement operators
    measOps[ct] = diagm(calData[:])

    #The data to invert

    tomoData = [tomoData;  real(datatemp[ct]["data"][1:end-nbrCalRepeats*(2^nbrQubits)])]

    #variance
    varData = [varData; datatemp[ct]["realvar"][1:end-nbrCalRepeats*(2^nbrQubits)]]

  end
  #Map each experiment to the appropriate readout pulse
  measOpMap = reshape(repmat(transpose(1:numMeas), nbrPulses^nbrQubits, 1), numMeas*(nbrPulses^nbrQubits), 1)

  measPulseMap = repmat((1:nbrPulses^nbrQubits), numMeas, 1)
  #Use a helper to get the measurement unitaries.
  measPulseUs = tomo_gate_set(nbrQubits, nbrPulses)
  varMat = diagm(varData[:])

  #Now call the inversion routines

  #First least squares
  rhoLSQ = QST_LSQ(tomoData, varMat, measPulseMap, measOpMap, measPulseUs, measOps, nbrQubits);
  #plotting to be implemented    pauliSetPlot(rho2pauli(rhoLSQ), newplot);

  #Now constrained SDP
  rhoSDP = QST_SDP(tomoData, varMat, measPulseMap, measOpMap, measPulseUs, measOps, nbrQubits);

  return (rhoLSQ, rhoSDP)
end

"""
    rho2pauli(rho)

Convert a density matrix to a Pauli set
"""
function rho2pauli(rho)
  nbrQubits = Int(log2(size(rho,1)));
  kronMat = zeros(Int, 4^nbrQubits, nbrQubits);
  for qubitct = 1:nbrQubits
    kronMat[:,qubitct] = reshape(repmat((1:4)', 4^(nbrQubits-qubitct), 4^(qubitct-1)), 4^nbrQubits, 1);
  end

  pauliVec = ones(4^nbrQubits,1);
    for ii = 1:4^nbrQubits
        pauliVec[ii] = real(trace(rho*allpaulis(1)[kronMat[ii,:][:]][1]));
    end

    pauliVec =  pauliVec[sortperm(map(weight,allpaulis(nbrQubits)))]
    pauliStrs = map(AbstractString, allpaulis(nbrQubits)[sortperm(map(weight,allpaulis(nbrQubits)))]) #reorder
    return pauliVec, pauliStrs
end
