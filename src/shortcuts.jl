using HDF5, Compat, KernelDensity

"""
    cal_data(data::Array; bit = 1, nqubits = 1, num_repeats = 2, rm_cals=true)

Normalize data using calibration points.  Note the data should not be complex.

# Arguments
  - `data::Array{Float64}`: data to normalize with cals at the end
  - `bit::Int`: 1, 2, 3... from LSB to MSB
  - `nqubits::Int`: number of qubits
  - `num_repeats::Int`: number of calibration points per computational state
  - `rm_cals::boolean`: remove calibration points from the end of the sequences

# Examples
```julia-repl
julia> data = vcat(rand(10,1), [0,0,1,1])
julia> Qlab.cal_data(data)
-0.037458  -0.763034  -0.830473  -1.69414  …  -1.25552  -2.45348  -2.45213
```
"""
function cal_data(data::Array;
                  bit = 1,
                  nqubits = 1,
                  num_repeats = 2,
                  rm_cals = true)
	zero_cal = mean(data[end-num_repeats*(2^nqubits)+1:end-num_repeats*(2^nqubits-1)])
	one_cal = mean(data[end-num_repeats*(2^nqubits-2^(bit-1))+1:end-num_repeats*(2^nqubits-2^(bit-1)-1)])
	scale_factor = -(one_cal - zero_cal)/2;
	if rm_cals
		data = data[1:end-(2^nqubits)*num_repeats]
	end
	data = (data .- zero_cal)/scale_factor .+ 1
end

"""
  cal_data(data; qubit="", cal0="0", cal1="1")

Normalize data with reference measurements defined in metadata. Format used in
Auspex.  Assume data with structure: [data, data, ...., cal0, cal0, ..., cal1,
cal1, cal1, ..., data, ..., cal0, ...].  The number of data and cals can be
different for each set, as long as they are contiguous.

# Arguments
  - `data`: dictionary (Auspex format)
  - `qubit::string`: qubit name
  - `cal0/1::string`: reference measurement for qubit in 0/1
"""

function cal_data(data::Tuple{Dict{String,Array{Any,N} where N},Dict{String,Any}};
                  qubit::String = "",
                  cal0::String = "0",
                  cal1::String = "1",
                  quad = :real)
    if length(collect(keys(data[1]))) == 1
        qubit = collect(keys(data[1]))[1]
    elseif isempty(qubit)
        error("More than one qubit. Specify qubit name and calibration labels")
    elseif ~(qubit in keys(data[1]))
        error("Qubit not found")
    end

    metadata = data[2][qubit]["meta_data"]
    if(length(metadata))>1
       error("Invalid metadata")
    end

    metadata_values = data[2][qubit]["meta_data"][(collect(keys(metadata)))[1]]

    data_out = []
    ind0 = findall(x -> x==cal0, metadata_values)
    ind1 = findall(x -> x==cal1, metadata_values)
    ind_data = findall(x -> x==metadata_values[1], metadata_values)
    ind0_edge = push!(filter(x -> ind0[x] != ind0[x+1]-1, 1:length(ind0)-1), length(ind0)) #find consecutive cals
    ind1_edge = push!(filter(x -> ind1[x] != ind1[x+1]-1, 1:length(ind1)-1), length(ind1))
    ind_data_edge = push!(filter(x -> ind_data[x] != ind_data[x+1]-1, 1:length(ind_data)-1), length(ind_data)) #find consecutive data. Assume that every sweep starts with data
    for s in 1:length(ind0_edge)
        ind0_s =  s>1 ? ind0[ind0_edge[s-1]+1:ind0_edge[s]] : ind0[1:ind0_edge[1]]
        ind1_s =  s>1 ? ind1[ind1_edge[s-1]+1:ind1_edge[s]] : ind1[1:ind1_edge[1]]
        ind_data_s =  s>1 ? ind_data[ind_data_edge[s-1]+1:ind_data_edge[s]] : ind_data[1:ind_data_edge[1]]
        data_s = eval(Symbol(quad)).(data[1][qubit][ind_data_s])
        zero_cal = mean(eval(Symbol(quad)).(data[1][qubit][ind0_s]))
        one_cal = mean(eval(Symbol(quad)).(data[1][qubit][ind1_s]))
        scale_factor = -(one_cal - zero_cal)/2
        data_s = (data_s .- zero_cal)/scale_factor .+ 1
        push!(data_out, data_s)
    end
      return data_out
end

"""
    get_fidelity(shots_0, shots_1; nbins=51, showPlot=false)

Get readout fidelity from single-shot measurements
"""
function get_fidelity(shots_0, shots_1, nbins = 51, showPlot = false)
  bins = range(min(minimum(shots_0),minimum(shots_1)), stop=max(maximum(shots_0),length=maximum(shots_1)),nbins)
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
  fidelity = 1-(1-maximum(abs.(cdf_0-cdf_1)))/2
  threshold = bins[argmax(abs.(cdf_0-cdf_1))]
  return fidelity, threshold
end

"""
    nanmean(vec)

Calculate the mean of vec, ignoring nan values.
"""
function nanmean(vec)
    return sum(.~isnan.(vec).*vec)/length(findall(.~isnan.(vec)))
end

"""
    get_extr_loc(M, dim; getmax=true)

Return index and value of max/min (getmax = True/False) in each column/row
(dim=1/2) of M
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
    get_feedback_fidelity(data::Int,
                          nqubits::Int,
                          nrounds::Int,
                          bit;
                          cal_repeats=2,
                          target_state=0)

Get fidelity for resetting 'bit' in a 'nqubits' register. Note that this
estimate assumes no measurement crosstalk.


# Arguments
  - `nrounds::Int`: number of feedback rounds
  - `bit`: 1, 2, 3... from LSB to MSB
  - `cal_repeats`: number of calibration points per computational state
  - `target_state`: 0/1 to target either qubit state
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

"""Unwrap angles by changing deltas between values to 2π complement.

Unwrap radian phase ϕ by changing absolute jumps greater than `discont` to
their 2π complement.

Args:
  ϕ: Input array.
  discont: Maximum discontinuity between values. Optional, default is π.

Returns:
  out: Output array of unwrapped phases.

See also:
  numpy.unwrap()
"""
function unwrap!(ϕ::Array{T}; discont=π) where T <: AbstractFloat
  Δ = diff(ϕ)
  Δmod = mod.(Δ .+ π, 2 * π) .- π
  Δmod[(Δmod .== -π) .& (Δ .> 0)] .= π
  ϕcorr = Δmod .- Δ
  ϕcorr[abs.(Δ) .< discont] .= 0
  return ϕ .+ vcat(0, cumsum(ϕcorr))
end

""" In-place version of unwrap. """
function unwrap(ϕ::Array{T}; discont=π) where T <: AbstractFloat
  return unwrap!(copy(ϕ), discont=discont)
end

function CR_hamiltonian(datapath,
                        filelist,
                        subdir=Dates.format(Dates.today(),"yymmdd"),
                        params_estimate = [5,5,0,0,0,1])
    IX = zeros(length(filelist))
    IY = zeros(length(filelist))
    IZ = zeros(length(filelist))
    ZX = zeros(length(filelist))
    ZY = zeros(length(filelist))
    ZZ = zeros(length(filelist))
    fit_curves = Dict()

    for (ct, filenum) in enumerate(filelist)
        data = load_data(datapath, filenum, datestr);

        norm_data = Qlab.cal_data(data[1], qubit="main")[1]
        data_len = Int(length(norm_data)/6)
        xvec = reshape(norm_data[1:3:end], data_len, 2)
        yvec = reshape(norm_data[2:3:end], data_len, 2)
        zvec = reshape(norm_data[3:3:end], data_len, 2);
        rvec = sqrt.(sum(xvec,2).^2+sum(yvec,2).^2+sum(zvec,2).^2)/2 # missing factor of 2 in PRA
        xpoints = data[2]["q1-main"][1]["points"][1:3:end-4][1:data_len]

        if ct>1
            params_estimate = tau_vec[1,:]
        end
        tau_vec = zeros(2,6)
        for ind=1:2
            result = optimize(x-> calculate_residueCR(x, xpoints, [xvec[:,ind],yvec[:,ind],zvec[:,ind]])[1], params_estimate)
            τ = Optim.minimizer(result)
            tau_vec[ind, :] = τ
            fit_curves[ind-1] = calculate_residueCR(τ, xpoints,  [xvec[:,ind],yvec[:,ind],zvec[:,ind]])[2]
        end
        IX[ct] = (tau_vec[1,1]+tau_vec[2,1])/2/(2*pi)
        IY[ct] = (tau_vec[1,2]+tau_vec[2,2])/2/(2*pi)
        IZ[ct] = (tau_vec[1,3]+tau_vec[2,3])/2/(2*pi)
        ZX[ct] = (tau_vec[1,1]-tau_vec[2,1])/2/(2*pi)
        ZY[ct] = (tau_vec[1,2]-tau_vec[2,2])/2/(2*pi)
        ZZ[ct] = (tau_vec[1,3]-tau_vec[2,3])/2/(2*pi)
        rfit = sqrt.(sum(hcat(fit_curves[0]+fit_curves[1]...).^2,2))/2
    end
end
