#collection of commonly used plots
function plot_ss_hists(shots_0, shots_1)
  #single-shot readout histograms
  hist_0 = kde(shots_0[:])
  hist_1 = kde(shots_1[:])
  fidelity = get_fidelity(shots_0, shots_1)
  sns.set(style="ticks")
  sns.set_style(Dict("xtick.direction" => "in"))
  sns.set_style(Dict("ytick.direction" => "in"))
  w=4
  figure(figsize=(w,w/1.4))
  plot(hist_0,label=L"$|0\rangle$")
  plot(hist_1,label=L"$|1\rangle$")
  ylabel("Counts")
  xlabel("Homodyne voltage (a.u.)")
  lgd = legend()
  annotate(@sprintf("Fid. = %0.2f", fidelity),
  xy=[0.1;0.7],
  xycoords="axes fraction",
  xytext=[0,10],
  textcoords="offset points",
  fontsize=10.0,
  ha="left",
  va="center")
  return hist_0, hist_1
end


function plot2D(data, quad = "real"; normalize=false)
  fig = figure("pyplot_surfaceplot",figsize=(5,3))
  ax = gca()
  ax[:ticklabel_format](useOffset=false)
  if quad == "real"
    data_quad = real(data["data"])
  elseif quad == "imag"
    data_quad = imag(data["data"])
  elseif quad == "amp"
    data_quad = abs(data["data"])
  end
  if normalize
    data_quad./=data_quad[:,1]
  end
  xpoints = repmat(data["xpoints"],1,length(data["ypoints"]))
  ypoints = repmat(data["ypoints"]',length(data["xpoints"]),1)
  pcolormesh(xpoints, ypoints, data_quad',cmap = "terrain")
  colorbar()
  xlabel(data["xlabel"])
  ylabel(data["ylabel"])
  xlim([minimum(xpoints),maximum(xpoints)])
  ylim([minimum(ypoints),maximum(ypoints)])
end
