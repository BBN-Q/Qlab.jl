using HDF5, Compat, KernelDensity

function cal_data(data; caltype="test", num_repeats = 2)
  zeroCal = mean(data[end-2*num_repeats+1:end-num_repeats]);
  piCal = mean(data[end-num_repeats+1:end]);

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
