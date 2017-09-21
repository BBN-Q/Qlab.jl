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

function plot1D(data, group = "main"; quad = "real", label_y = "V (a.u.)", normalize = false, bit = 1, nqubits = 1, num_repeats = 2, save_fig = "png")
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
   title(get_partial_filename(data[3]["filename"]))
   ~isempty(save_fig) && savefig(string(splitext(data[3]["filename"])[1],'-',group,'.', save_fig))

  return xpoints_values, data_values
end

function plot2D(data, group = "main"; quad = "real", transpose = false, normalize = false, vmin = NaN, vmax = NaN, cmap = "terrain", show_plot = true, save_fig = "png")
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
    title(get_partial_filename(data[3]["filename"]))
    ~isempty(save_fig) && savefig(string(splitext(data[3]["filename"])[1],'-',group,'.', save_fig))
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

function plot_multi(data, group = "main"; quad = "real", offset = 0.0, cals = false, show_legend = true,
  cal0::String = "0", cal1::String = "1", fit_name = "", fit_param_name = "T", save_fig = "png")
  data_values = data[1][group]["Data"]
  xpoints = data[2][group][2]
  ypoints = data[2][group][1]
  xpoints_values = xpoints["points"]
  ypoints_values = ypoints["points"]
  if isempty(fit_name)
    figure(figsize= (3.5,3))
  else
    figure(figsize = (8,3))
    subplot(1,2,1)
    Tvec = zeros(length(ypoints_values))
    dTvec = zeros(length(ypoints_values))
    fit_function = eval(parse(string("fit_", fit_name)))
  end
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
  ax = gca()
  for k in 1:length(data_quad)
    xpts = xpoints_values[1:length(data_quad[k])]
    plot(xpts, data_quad[k] + offset*k, label=string(ypoints_values[k]), linewidth = convert(Int,isempty(fit_name)), marker = "o", markersize = 4)
    if ~isempty(fit_name)
      fit_result = (fit_function)(xpts, data_quad[k])
      k==1 && println("Fitting to model: ", fit_result.model_str)
      plot(xpts, fit_result.fit_curve(xpts),label="fit", color=ax[:lines][end][:get_color](), linewidth=1)
      Tvec[k] = fit_result.fit_params[fit_param_name]
      dTvec[k] = fit_result.errors[fit_param_name]
    end
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
    plr = subplot(1,2,2)
    errorbar(ypoints_values, Tvec, dTvec, marker = "o", markersize=4)
    ylabel(fit_param_name) #TODO: unit? Proper parameter name, e.g., T1?
    xlabel(ypoints["name"])
    plr[:axis](ymin = 0)
    subplots_adjust(wspace=0.3)
    return xpoints_values[1:length(data_quad[1])], data_quad, (Tvec, dTvec)
  end
  title(get_partial_filename(data[3]["filename"]))
  ~isempty(save_fig) && savefig(string(splitext(data[3]["filename"])[1],'-',group,'.', save_fig))
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

function get_partial_filename(filename, num_dirs = 2)
  cur_path = ""
  for n=1:num_dirs+1
    filename, filename_dir = splitdir(filename)
    cur_path = n==1? filename_dir : joinpath(filename_dir, cur_path)
  end
  return cur_path
end
