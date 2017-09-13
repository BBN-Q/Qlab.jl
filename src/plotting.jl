using PyPlot, KernelDensity, Formatting, Seaborn
# using Formatting
#collection of commonly used plots
function plot_ss_hists(shots_0, shots_1)
  #single-shot readout histograms
  hist_0 = kde(shots_0[:])
  hist_1 = kde(shots_1[:])
  fidelity, _ = get_fidelity(shots_0, shots_1)
  w=4
  figure(figsize=(w,w/1.4))
  plot(hist_0.x, hist_0.density/sum(hist_0.density),label=L"$|0\rangle$")
  plot(hist_1.x, hist_1.density/sum(hist_1.density),label=L"$|1\rangle$")
  ylabel("Fraction of counts")
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

function plot1D(data, group = "main"; quad = "real", label_y = "V (a.u.)", normalize = false, bit = 1, nqubits = 1, num_repeats = 2)
  fig = figure(figsize=(3,3))
  data_values = data[1][group]["Data"]
  xpoints = data[2][group][1]
  xpoints_values = xpoints["points"]
  quad_f = Dict("real"=> real, "imag"=> imag, "amp"=> abs, "abs" => abs)
  data_quad = quad_f[quad].(data_values)
  if normalize
    data_quad = Qlab.cal_data(data_quad)
    xpoints_values = xpoints_values[1:end-num_repeats*(2^nqubits)]
    label_y = L"\langle Z\rangle"
end
   plot(xpoints_values, data_quad)
   label_x = xpoints["name"]
   if xpoints["unit"] != "None"
    label_x = string(label_x, " (", xpoints["unit"], ")" )
   end
   xlabel(label_x)
   ylabel(label_y)
  return xpoints_values, data_values
end

function plot2D(data, group = "main"; quad = "real", transpose = false, normalize = false, vmin = NaN, vmax = NaN, cmap = "terrain", show_plot = true)
  data_values = data[1][group]["Data"]
  xpoints = data[2][group][1]
  ypoints = data[2][group][2]
  xpoints_values = xpoints["points"]
  ypoints_values = ypoints["points"]
  quad_f = Dict("real"=> real, "imag"=> imag, "amp"=> abs, "abs" => abs)
  data_quad = quad_f[quad].(data_values)
  xpoints_grid = repmat(xpoints_values', length(ypoints_values), 1)
  ypoints_grid = repmat(ypoints_values, 1, length(xpoints_values))
  data_grid = reshape(data_quad, length(ypoints_values), length(xpoints_values))
  if normalize == 1
   data_grid = (data_grid'./data_grid[1,:])'
  elseif normalize == 2
   data_grid./=data_grid[:,1]
  end
  label_x = xpoints["name"]
  contains(label_x, "_metadata") && (label_x = split(label_x, "_metadata")[1])
  if xpoints["unit"] != "None"
    label_x = string(label_x, " (", xpoints["unit"], ")" )
  end
  label_y = ypoints["name"]
  contains(label_y, "_metadata") && (label_y = split(label_y, "_metadata")[1])
  if ypoints["unit"] != "None"
    label_y = string(label_y, " (", ypoints["unit"], ")" )
  end
  if isnan(vmin)
    vmin = minimum(data_grid)
  end
  if isnan(vmax)
    vmax = maximum(data_grid)
  end
  if show_plot
    fig = figure("pyplot_surfaceplot",figsize=(3,3))
    ax = gca()
    ax[:ticklabel_format](useOffset=false)
    if transpose
      pcolormesh(ypoints_grid, xpoints_grid, data_grid, cmap = cmap, vmin = vmin, vmax = vmax)
      ylabel(label_x)
      xlabel(label_y)
    else
      pcolormesh(xpoints_grid, ypoints_grid, data_grid, cmap = cmap, vmin = vmin, vmax = vmax)
      xlabel(label_x)
      ylabel(label_y)
    end
    colorbar()
  end
  return xpoints_values, ypoints_values, data_grid
end

function reshape2D(data, group = "main"; quad = "real", normalize = false)
  return plot2D(data, group; quad = quad, normalize = normalize, show_plot = false)
end


"""
  plot_multi

Plot multiple 1D traces on top of each other.
group: data group name
Optional arguments
-------------------------
quad: quadrature: real, imag, or abs
offset: vertical offset between curves
cals: normalize to 0/1 using metadata
show_legend: show legend in plot
"""

function plot_multi(data, group = "main"; quad = "real", offset = 0.0, cals = false, show_legend = true, cal0::String = "0", cal1::String = "1", fit_name = "")
  figure(figsize = (isempty(fit_name)? 6 : 3,3))
  subplot(1, isempty(fit_name)? 1 : 2, 1)
  data_values = data[1][group]["Data"]
  xpoints = data[2][group][2]
  ypoints = data[2][group][1]
  xpoints_values = xpoints["points"]
  ypoints_values = ypoints["points"]
  quad_f = Dict("real"=> real, "imag"=> imag, "amp"=> abs, "abs" => abs)
  data_quad = quad_f[quad].(data_values)
  if cals
    data_quad = cal_data(data[1], qubit=group)
    data_quad = quad_f[quad].(data_quad)
  else
    # reshape to array of array
    data_quad = reshape(data_quad, length(xpoints_values), length(ypoints_values))
    data_quad = [data_quad[:,k] for k in 1:size(data_quad,2)]
  end
  for k in 1:length(data_quad)
    plot(xpoints_values[1:length(data_quad[k])], data_quad[k] + offset*k, label=string(ypoints_values[k]))
  end
  label_x = xpoints["name"]
  contains(label_x, "_metadata") && (label_x = split(label_x, "_metadata")[1])
  if xpoints["unit"] != "None"
    label_x = string(label_x, " (", xpoints["unit"], ")" )
  end
  xlabel(label_x)
  if cals
    ylabel("<Z>")
  else
    ylabel("Voltage (a.u.)")
  end
  show_legend && (legend())

  if ~isempty(fit_name)
    subplot(1,2,2)
    plot_fit_multi(xpoints_values[1:length(data_quad[1])], data_quad, fit_name)
  end
  return xpoints_values[1:length(data_quad[1])], data_quad
end

function plot2D_matlab(data, quad = "real"; normalize=false)
  fig = figure("pyplot_surfaceplot",figsize=(5,3))
  ax = gca()
  ax[:ticklabel_format](useOffset=false)
  if quad == "real"
    data_quad = real(data["data"])
  elseif quad == "imag"
    data_quad = imag(data["data"])
  elseif quad == "amp"
    data_quad = abs.(data["data"])
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

function pauli_set_plot(rho; rho_ideal=[], fig_width=5, fig_height=3.5, bar_width=0.6)
    pauli_vec, pauli_ops = rho2pauli(rho)
    figure(figsize=(fig_width,fig_height))
    ind = 1:length(pauli_vec)
    if ~isempty(rho_ideal)
        pauli_vec_ideal, _  = rho2pauli(rho_ideal)
        bar(ind, pauli_vec_ideal, bar_width, color="green", label=L"$\rho_{ideal}$")
    end
    bar(ind, pauli_vec, bar_width, label=L"$\rho_{actual}$")
    xticks(ind + bar_width/2., map(string,pauli_ops))
    ylim([-1.05,1.05])
    xlabel("Pauli operator")
    ylabel("Expectation value")
    title("State tomography")
    if ~isempty(rho_ideal)
        legend()
    end
end


function annotate_plot(message, vals...; coords = [0.75, 0.9], fontsize = 10.0)
  annotate(format(message, vals...),
    xy= coords,
    xycoords="axes fraction",
    xytext=[10,15],
    textcoords="offset points",
    fontsize=fontsize,
ha="left",
va="center")
end

function plot_fit_multi(xpts, ypts, function_name, iter = [])
  Tvec = zeros(length(ypts))
  dTvec = zeros(length(ypts))
  fit_function = eval(parse(string("fit", function_name)))
  figure(figsize = (8,3))
  subplot(1,2,1)
  ax = gca()
  for k=1:length(ypts)
      fit_result = (fit_function)(xpts,ypts[k])
      plot(xpts, ypts[k], label="data", marker = "o", linewidth = 0, markersize=2)
      plot(xpts, fit_result.fit_curve(xpts),label="fit", color=ax[:lines][end][:get_color](), linewidth=1)
      Tvec[k] = fit_result.fit_params["T"]
      dTvec[k] = fit_result.errors["T"]
      #Qlab.annotate_plot(L"$T_1 = {1:.1f}\pm{2:.1f}\,\mu s$", fit_result.fit_params["T"], fit_result.errors["T"], coords = (0.4, 0.7))
      xlabel(L"Time ($\mu$s)")
      ylabel(L"\langle Z \rangle")
  end
  plr = subplot(1,2,2)
  if isempty(iter)
    iter = 1:length(ypts)
    xlabel_name = "Repeat"
  end
  isempty(iter) && (iter = 1:length(ypts))
  errorbar(iter, T1vec, dT1vec, marker = "o", markersize=4)
  ylabel(L"$T_1$ ($\mu$ s)")
  xlabel(xlabel_name)
  plr[:axis](ymin = 0)
end
