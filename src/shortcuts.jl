using HDF5, Compat, KernelDensity

function cal_data(data; caltype="test", num_repeats = 2)

  if caltype == "a"
    zeroCal = (mean(data[end-2*num_repeats+1:4:end])+mean(data[end-2*num_repeats+2:4:end]))/2;
    piCal = (mean(data[end-2*num_repeats+3:4:end])+mean(data[end-2*num_repeats+4:4:end]))/2;
  elseif caltype == "b"
    zeroCal = mean(data[end-2*num_repeats+1:end-num_repeats]);
    piCal = mean(data[end-num_repeats+1:end]);
  else
    zeroCal = mean(data[end-2*num_repeats+1:end-num_repeats]);
    piCal = mean(data[end-num_repeats+1:end]);
  end

  scaleFactor = -(piCal - zeroCal)/2;
  data = data[1:end-2*num_repeats];
  data = (data - zeroCal)/scaleFactor + 1;
end

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
  return fidelity
end

function nanmean(vec)
    return sum(~isnan(vec).*vec)/length(find(~isnan(vec)))
end

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

function reshape_cells(data)
  measnum=length(data);
  cellsize=size(data[1]["data"])
  repeatnum =cellsize[1];
  alldata = Dict()
  for k=1:measnum
          alldata[k]=zeros(1,0);
      for r=1:repeatnum
              tempdata = data[k]["data"]
              alldata[k]=hcat(alldata[k], tempdata[r,:]);
      end
  end
  return alldata
end

function get_max_loc(M, dim)
  #return index and value of max in each column (1)/row (2)
  _, aindx = findmax(M, dim);

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

function get3pops(data)
#get transmon populations after an experiment lternating Id, Pi_01, Pi_12 as tomography pulses
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
