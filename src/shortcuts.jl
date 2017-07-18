using HDF5, Compat, KernelDensity

"""
    cal_data(data; bit, nqubits, num_repeats)

Normalize data using calibration points
bit: 1, 2, 3... from LSB to MSB
nqubits: number of qubits
num_repeats: number of calibration points per computational state
"""
function cal_data(data::Array; bit = 1, nqubits = 1, num_repeats = 2)
	zero_cal = mean(data[end-num_repeats*(2^nqubits)+1:end-num_repeats*(2^nqubits-1)])
	one_cal = mean(data[end-num_repeats*(2^nqubits-2^(bit-1)):end-num_repeats*(2^nqubits-2^(bit-1)-1)-1])
	scale_factor = -(one_cal - zero_cal)/2;
	data = data[1:end-2*num_repeats]
	data = (data - zero_cal)/scale_factor + 1
end

"""
  cal_data(data; qubit, cal0, cal1)

Normalize data with reference measurements defined in metadata. Format used in Auspex.

data: dictionary (Auspex format)
qubit: qubit name
cal0/1: reference measurement for qubit in 0/1
"""

function cal_data(data::Dict{String,Dict{String,Array{Any,1}}}; qubit::String = "", cal0::String = "0", cal1::String = "1")
  if length(collect(keys(data))) == 1
    qubit = collect(keys(data))[1]
  elseif isempty(qubit)
    error("More than one qubit. Specify qubit name and calibration labels")
  elseif ~(qubit in keys(data))
    error("Qubit not found")
  end

  metadata_key = []
  for k in keys(data[qubit])
    if contains(k, "metadata")
      push!(metadata_key, k)
    end
  end
  if length(metadata_key) != 1
      error("Invalid metadata")
  end
  ind0 = find(x -> x==parse(UInt8, cal0, 2), data[qubit][metadata_key[1]])
  ind1 = find(x -> x==parse(UInt8, cal1, 2), data[qubit][metadata_key[1]])
  zero_cal = mean(data[qubit]["Data"][ind0])
  one_cal = mean(data[qubit]["Data"][ind1])
  scale_factor = -(one_cal - zero_cal)/2
  data = data[qubit]["Data"][find(x -> x == data[qubit][metadata_key[1]][1], data[qubit][metadata_key[1]])] #assumes calibrations at the end only
  data = (data - zero_cal)/scale_factor + 1
end

"""
    get_fidelity(shots_0, shots_1; nbins, showPlot)

Get readout fidelity from single-shot measurements
"""
function get_fidelity(shots_0, shots_1, nbins = 51, showPlot = false)
  bins = linspace(min(minimum(shots_0),minimum(shots_1)),max(maximum(shots_0),maximum(shots_1)),nbins)
  hist_0 = kde(shots_0[:],bins)
  hist_1 = kde(shots_1[:],bins)
  cdf_0 = cumsum(hist_0.density)/cumsum(hist_0.density)[end]
  cdf_1 = cumsum(hist_1.density)/cumsum(hist_1.density)[end]
  if showPlot
    plot(hist_0.x, cdf_0)
    plot(hist_1.x, cdf_1)
    ylabel("Cumulative counts")
    xlabel("Homodyne voltage (a.u.)")
  end
  fidelity = 1-(1-maximum(abs(cdf_0-cdf_1)))/2
  threshold = bins[indmax(abs(cdf_0-cdf_1))]
  return fidelity, threshold
end

"""
    nanmean(vec)

Calculate the mean ignoring nan values
"""
function nanmean(vec)
    return sum(~isnan(vec).*vec)/length(find(~isnan(vec)))
end

"""
    get_extr_loc(M, dim; getmax)

Return index and value of max/min (getmax = True/False) in each column/row (dim=1/2) of M
"""
function get_extr_loc(M, dim, getmax = true)
  #return index and value of max/min in each column (1)/row (2)
  if getmax
    _, aindx = findmax(M, dim)
  else
    _, aindx = findmin(M, dim)
  end
  msize=size(M);
  col = zeros(msize[-dim+3], 1);
  for i=1:msize[-dim+3]
    if dim == 1
      col[i], _ = ind2sub(msize, aindx[i])
    elseif dim == 2
      _, col[i] = ind2sub(msize, aindx[i])
    end
  end
  return convert(Array{Int8,2}, col), M[aindx]
end

"""
    get3pops(data)

get transmon populations after an experiment alternating Id, pi_01, pi_12 as tomography pulses
"""
function get3pops(data)
  calRepeat = 2;
  pvec3 = zeros(3,Int((length(data)-3*calRepeat)/3));
  pvec2 = zeros(3,Int((length(data)-3*calRepeat)/3));
  pvec2temp = zeros(2,Int((length(data)-3*calRepeat)/3));
  v0 = mean(data[end-3*calRepeat+1:end-2*calRepeat]);
  v1 = mean(data[end-2*calRepeat+1:end-calRepeat]);
  v2 = mean(data[end-calRepeat+1:end]);
  m0 = data[1:3:end]
  m1 = data[2:3:end]
  m2 = data[3:3:end]

  M = [v0 v1 v2; v1 v0 v2; v0 v2 v1];
  for c = 1:Int((length(data)-3*calRepeat)/3)
    pvec3[:,c] = M\[m0[c];m1[c];m2[c]]
  end

  #fixed p0+p1+p2 = 1
  M2 = [v0-v2 v1-v2; v1-v2 v0-v2];
  for c = 1:Int((length(data)-3*calRepeat)/3)
    pvec2temp[:,c] = M2\[m0[c]-v2;m1[c]-v2];
  end
  pvec2 = [pvec2temp; 1-sum(pvec2temp,1)];
  return (pvec3, pvec2)
end

"""
    get_feedback_fidelity(data, nqubits, nrounds, bit; cal_repeats = 2, target_state = 0)

get fidelity for resetting 'bit' in a 'nqubits' register. Note that this estimate assumes no measurement crosstalk.
nrounds: feedback rounds
bit: 1, 2, 3... from LSB to MSB
cal_repeats: number of calibration points per computational state
target_state: 0/1 to target either qubit state
"""
function get_feedback_fidelity(data, nqubits, nrounds, bit; cal_repeat = 2, target_state = 0)
    zero_cal = mean(data[end-cal_repeat*(2^nqubits)+1:end-cal_repeat*(2^nqubits-1)])
    one_cal = mean(data[end-cal_repeat*(2^nqubits-2^(nqubits-bit))+1:end-cal_repeat*(2^nqubits-1-2^(nqubits-bit))])
    scale_factor = -(one_cal - zero_cal)/2;
    data = data[1:end-cal_repeat*2^nqubits]
    normdata = (data - zero_cal)./scale_factor + 1
    fidmat = zeros(2^nqubits, nrounds)
    for state_ind = 1:size(fidmat,1)
        for round_ind = 1:size(fidmat,2)
            fidmat[state_ind, round_ind] = (1+(-1)^target_state*normdata[(state_ind-1)*(nrounds+1)+round_ind+1])/2
        end
    end
    return fidmat
end
