using StatsBase

"""
   digitize(data, cal0, cal1, mode=:equal)

  `data`  --  iterable holding floating point values to be digitized

  `cal0`  --  vectors of floating point values corresponding to samples of
  `cal1`      voltages from measurements of the 0 and 1 states

  `mode`  --  mode of operation for thresholding. Supported options are `:equal` and `:max`

  Translate floating point values in `data` to binary outcomes based on
calibration samples for 0 and 1 outcomes. It is assumed that the number of samples
in `cal0` and `cal1` is large enough that the lower and upper quartile estimate are
accurate. The output is a `(ddata, pe0, pe1)`, where

  `ddata` -- digitized version of ddata (0s and 1s as `Int`s)

  `pe0`   -- probability of error when preparing a 0

  `pe1`   -- probability of error when preparing a 1

  `thr`   -- chosen threshold

  The `:equal` mode of operation tries to equalize the conditional error probabilities for the
two possible outcomes, and it is generally less sensitive to statistical fluctuations. The
`:max` mode of operation tries to minimize the sum of error probabilities, but is is much more
sensitive to statistical fluctuations, and is therefore not recommended unless the ensemble
size for the calibrations is very large.

"""
function digitize(data, cal0::Vector{Float64}, cal1::Vector{Float64}; mode = :equal, thr = NaN)
    minA, maxA = extrema([quantile(cal0,[.25,.75]); quantile(cal1,[.25,.75])])
    Arange = linspace(minA,maxA,1000)
    if isnan(thr)
      thr = if mode==:max
              Arange[indmax(abs(ecdf(cal0)(Arange)-ecdf(cal1)(Arange)))]
            elseif mode==:equal
              Arange[indmin(abs(ecdf(cal0)(Arange)-1+ecdf(cal1)(Arange)))]
            else
              error("Unsupported digitization mode")
            end
    end

    ddata = data .> thr
    dcal0 = cal0 .> thr
    dcal1 = cal1 .> thr

    if sum(dcal0)/length(cal0) > 1/2
        ddata = ~ddata
        dcal0 = ~dcal0
        dcal1 = ~dcal1
    end

    return round(Int,ddata), sum(dcal0)/length(cal0), 1.0-sum(dcal1)/length(cal1), thr
end

"""
    digitize_2D(k00, k01, k10, k11, n)

Digitize readout voltages in a n x n grid, using the estimated
distributions for the computational states k00, k01, 01, k11
"""
function digitize_2D(k00,k01,k10,k11,n)
  xmin = minimum([k00.x; k01.x; k10.x; k11.x])
  xmax = maximum([k00.x; k01.x; k10.x; k11.x])
  ymin = minimum([k00.y; k01.y; k10.y; k11.y])
  ymax = maximum([k00.y; k01.y; k10.y; k11.y])
  xpoints = linspace(xmin,xmax,n)
  ypoints = linspace(ymin,ymax,n)
  assignmat = zeros(n,n)
  distrs = [k00 k01 k10 k11]
  for yi in 1:n
    y_ind = zeros(Int,4)
    for k=1:4
      y_ind[k] = max(1,searchsortedlast(distrs[k].y, ypoints[yi]))
    end
    for xi in 1:n
      x_ind = zeros(Int,4)
      for m=1:4
        x_ind[m] = max(1,searchsortedlast(distrs[m].x, xpoints[xi]))
      end
      highest_kde_ind =  indmax([distrs[1].density[x_ind[1],y_ind[1]] distrs[2].density[x_ind[2],y_ind[2]] distrs[3].density[x_ind[3],y_ind[3]] distrs[4].density[x_ind[4],y_ind[4]]])
      if distrs[highest_kde_ind].density[x_ind[highest_kde_ind],y_ind[highest_kde_ind]]>0
        assignmat[xi, yi] = highest_kde_ind
      end
    end
  end
  return xpoints, ypoints, assignmat
end

"""
  smooth_2D_map(digmat; dig_thr)

Smooth a 2D classification map by doing a majority vote over the nearest neighbors.
dig_thr is the minimum number of equal nearest neighbors to flip a state
Example: point digmat[i,j] is surrounded by 2 points equal to 1, 3 equal to 2, 3 equal to 3, 1 equal to 4.
If dig_thr <=3, digmat[i,j] will be set to a value chosen randomly between 2 and 3. dig_thr
is considered to be 1 if digmat[i,j] = 0 (unassigned)
"""
function smooth_2D_map(digmat; dig_thr = 8)
  roundvec = zeros(Int, 8) #8 nearest neighbors
  nn_count = zeros(Int, 4) #counts/state
  for n in 2:size(digmat,2)-1
    for m in 2:size(digmat,1)-1
      roundvec = [digmat[m-1,n] digmat[m+1,n] digmat[m, n-1] digmat[m, n+1] digmat[m-1, n+1] digmat[m-1, n-1] digmat[m+1, n+1] digmat[m+1, n-1]]
      for k = 1:4
        nn_count[k] = count(x -> x == k, roundvec)
      end
      highest_nn = maximum(nn_count)
      indmaxs = find(x -> x == highest_nn, nn_count) #find all the dominant colors
      indmax = indmaxs[rand(1:end)] #break ties randomly
      if digmat[m,n] == 0 || maximum(nn_count) >= dig_thr #flip states if surrounded
        #by dig_thr points of the same color, or if the state is unassigned
        digmat[m,n] = indmax
      end
    end
  end
  return digmat
end

"""
    get_fidelities_2D(distrs, xpoints, ypoints, digmat)

Calculate assignment fidelities in a 2-qubit space
distrs: distributions for the 4 computational states
digmat: digitization map with grid [xpoints; ypoints]
"""
function get_fidelities_2D(distrs, xpoints, ypoints, digmat)
  probmat = zeros(4,4)
  for k = 1:4
    #ensure that the distributions are normalized
    distrs[k].density/=sum(distrs[k].density)
    for (n, yi) in enumerate(distrs[k].y)
      #search in digmat
      y_ind = max(1,searchsortedlast(ypoints, yi))
      for (m, xi) in enumerate(distrs[k].x)
        x_ind = max(1,searchsortedlast(xpoints, xi))
        for d=1:4
          probmat[k,d] += distrs[k].density[m,n]*(digmat[x_ind,y_ind] == d)
        end
      end
    end
  end
  @printf("Assign fid. P_00 = %.2f\n", probmat[1,1])
  @printf("Assign fid. P_01 = %.2f\n", probmat[2,2])
  @printf("Assign fid. P_10 = %.2f\n", probmat[3,3])
  @printf("Assign fid. P_11 = %.2f\n\n", probmat[4,4])
  p_0x = sum(probmat[1:2,1:2])/2 #assignment prob. for A1 independent of A2
  p_1x = sum(probmat[3:4,3:4])/2
  p_x0 = (probmat[1,1]+probmat[1,3]+probmat[3,1]+probmat[3,3])/2
  p_x1 = (probmat[2,2]+probmat[2,4]+probmat[4,2]+probmat[4,4])/2
  @printf("Single-qubit fid. A1 = %.2f\n", (p_0x+p_1x)/2)
  @printf("Single-qubit fid. A2 = %.2f\n", (p_x0+p_x1)/2)
  return probmat
end

"""
    single_shot_fidelities_2D(filenum, datestr; ch_a1, ch_a2, grid_dim, nn_quorum, showPlot)

Wrapper to load and digitize 2-q single-shot data.
Data format: single-shot alternating between 00, 01, 10, 11 in channels {ch_a1, ch_a2}
grid_dim: number of bins / qubit
nn_quorum: number of votes required for the nearest neighbors to flip a value in the assignment map
"""
function single_shot_fidelities_2D(filenum, datestr; ch_a1 = 1, ch_a2 = 2, grid_dim=51, nn_quorum = 6, showPlot = false)
  data=load_data(datapath, filenum, datestr)[2];
  a1 = real(data[ch_a1]["data"])
  a2 = real(data[ch_a2]["data"])
  kde_00 = kde((a1[1:4:end],a2[1:4:end]))
  kde_01 = kde((a1[2:4:end],a2[2:4:end]))
  kde_10 = kde((a1[3:4:end],a2[3:4:end]))
  kde_11 = kde((a1[4:4:end],a2[4:4:end]))
  xpts,ypts,assign_mat = digitize_2D(kde_00,kde_01,kde_10,kde_11, grid_dim)
  smoothed_mat = smooth_2D_map(assign_mat; dig_thr = nn_quorum);
  fidelity_mat = get_fidelities_2D([kde_00 kde_01 kde_10 kde_11], xpts, ypts, smoothed_mat)
  if showPlot
    fig = figure("pyplot_surfaceplot",figsize=(5,4))
    xpoints = repmat(xpts,1,length(ypts))
    ypoints = repmat(ypts',length(xpts),1)
    ax = fig[:add_subplot](1,1,1)
    pcolormesh(xpts, ypts, smoothed_mat',cmap = PyPlot.get_cmap("Accent", 5),vmin=0,vmax=4)
    colorbar(ticks=0:4)
    xlabel(L"$A_1$ (V)")
    ylabel(L"$A_2$ (V)")
  end
  return smoothed_mat, fidelity_mat
end
