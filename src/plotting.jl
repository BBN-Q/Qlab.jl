using KernelDensity, Formatting, PyPlot, Seaborn, LinearAlgebra
# using Formatting
#collection of commonly used plots
function plot_ss_hists(shots_0, shots_1)
  #single-shot readout histograms
  hist_0 = kde(shots_0[:])
  hist_1 = kde(shots_1[:])
  fidelity, threshold = get_fidelity(shots_0, shots_1)
  println("Treshold = ", threshold)
  w=4
  figure(figsize=(w,w/1.4))
  plot(hist_0.x, hist_0.density/sum(hist_0.density),label=L"$|0\rangle$")
  plot(hist_1.x, hist_1.density/sum(hist_1.density),label=L"$|1\rangle$")
  ylabel("Fraction of counts")
  xlabel("Homodyne voltage (a.u.)")
  lgd = legend()
  annotate(sprint(show, "Fid. = $(fidelity)"),
  xy=[0.1;0.7],
  xycoords="axes fraction",
  xytext=[0,10],
  textcoords="offset points",
  fontsize=10.0,
  ha="left",
  va="center")
  return hist_0, hist_1
end

function plot1D(data, group = "main"; quad = :real, label_y = "V (a.u.)", cals = false, cal0::String = "0", cal1::String = "1", fit_name = "", save_fig = "png", doplot=true, fig = nothing)
  if fig == nothing && doplot == true
      fig = figure(figsize=(3,3))
  end
  xpoints = data[2][group]["axes"]
  xpoints_values = collect(values(xpoints))[1]
  if cals
    data_values = cal_data(data, qubit=group, quad = quad, cal0=cal0, cal1=cal1)[1]
    label_y = L"\langle Z\rangle"
  else
    data_values = dropdims(eval(quad).(data[1][group]),dims=1)
    label_y = string(quad, "(Voltage)")
  end
  xpoints_values = xpoints_values[1:length(data_values)]
  if doplot
      plot(xpoints_values, data_values, marker=".")
  end
  if ~isempty(fit_name)
    fit_function = eval(Meta.parse(string("fit_", fit_name)))
    fit_result = nothing
    try
        fit_result = (fit_function)(xpoints_values, data_values)
        println("Fitting to model: ", fit_result.model_str)
        ax = gca()
        if doplot
            plot(xpoints_values, fit_result.fit_curve(xpoints_values),label="fit", color=ax.lines[end].get_color(), linewidth=1)
        end
    catch
        println("Fit failed")
        return (xpoints_values, data_values, nothing)
    end
  end
  label_x = collect(keys(xpoints))[1]
  units = data[2][group]["units"][label_x]
  if units != nothing
    label_x = string(label_x, " (", units, ")" )
  end
  if doplot
      xlabel(label_x)
      ylabel(label_y)
      title(get_partial_filename(data[2][group]["filename"]))
      ~isempty(save_fig) && savefig(string(splitext(data[2][group]["filename"])[1],'-',group,'.', save_fig), bbox_inches = "tight")
  end
  if isempty(fit_name)
    return (xpoints_values, data_values)
  else
    return (xpoints_values, data_values, fit_result)
  end
end

function plot2D(data, group = "main"; quad = :real, transpose = false, normalize = false, vmin = NaN, vmax = NaN, cmap = "terrain", show_plot = true, save_fig = "png")
  data_values = data[1][group]
  axes = data[2][group]["axes"]
  xpoints_values = collect(values(axes))[1]
  ypoints_values = collect(values(axes))[2]
  data_quad = eval(quad).(data_values)
  if normalize == 1
   data_quad = (data_quad'./data_quad[1,:])'
  elseif normalize == 2
   data_quad./=data_quad[:,1]
  end
  label_x = collect(keys(axes))[1]
  occursin("_metadata", label_x) && (label_x = split(label_x, "_metadata")[1])
  units_x = data[2][group]["units"][label_x]
  if units_x != nothing
    label_x = string(label_x, " (", units_x, ")" )
  end
  label_y = collect(keys(axes))[2]
  units_y= data[2][group]["units"][label_y]
  occursin("_metadata", label_y) && (label_y = split(label_y, "_metadata")[1])
  if units_y != nothing
    label_y = string(label_y, " (", units_y, ")" )
  end
  if isnan(vmin)
    vmin = minimum(data_quad)
  end
  if isnan(vmax)
    vmax = maximum(data_quad)
  end
  if show_plot
    fig = figure("pyplot_surfaceplot",figsize=(4,4))
    ax = gca()
    ax[:ticklabel_format](useOffset=false)
    if transpose
      pcolormesh(ypoints_values, xpoints_values, data_quad', cmap = cmap, vmin = vmin, vmax = vmax)
      ylabel(label_x)
      xlabel(label_y)
    else
      pcolormesh(xpoints_values, ypoints_values, data_quad, cmap = cmap, vmin = vmin, vmax = vmax)
      xlabel(label_x)
      ylabel(label_y)
    end
    colorbar()
    if !haskey(data[2][group],"filename")
        data[2][group]["filename"] = "None"
    end
    title(get_partial_filename(data[2][group]["filename"]))
    ~isempty(save_fig) && savefig(string(splitext(data[2][group]["filename"])[1],'-',group,'.', save_fig), bbox_inches = "tight")
  end
  return xpoints_values, ypoints_values, data_quad
end

function reshape2D(data, group = "main"; quad = :real, normalize = false)
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
fit_name: name of fit function as in fit_###
fit_param_name: fit parameter of interest (output)
save_fig: format of the figure to be saved (empty string to disable)
"""

function plot_multi(data, group = "main"; quad = :real, offset = 0.0, cals = false, show_legend = true,
  cal0::String = "0", cal1::String = "1", fit_name = "", fit_param_name = "T", save_fig = "png")
  data_values = data[1][group]
  axes = data[2][group]["axes"]
  xpoints_values = collect(values(axes))[1]
  ypoints_values = collect(values(axes))[2]
  Tvec = zeros(length(ypoints_values))
  dTvec = zeros(length(ypoints_values))
  if isempty(fit_name)
    fig = figure(figsize= (3.5,3))
  else
    fig = figure(figsize = (8,3))
    subplot(1,2,1)
    fit_function = eval(Meta.parse(string("fit_", fit_name)))
  end
  if cals
    data_values = cal_data(data[1], qubit=group, quad=quad, cal0=cal0, cal1=cal1)
  else
    data_values = eval(quad).(data[1][group])
    # reshape to array of array
    data_values = reshape(data_values, length(xpoints_values), length(ypoints_values))
    data_values = [data_values[:,k] for k in 1:size(data_values,2)]
  end
  ax = gca()
  for k in 1:length(data_values)
    xpts = xpoints_values[1:length(data_values[k])]
    plot(xpts, data_values[k] .+ offset*k, label=string(round(ypoints_values[k],digits=3)), linewidth = convert(Int,isempty(fit_name)), marker = "o", markersize = 4)
    if ~isempty(fit_name)
      try
          fit_result = (fit_function)(xpts, data_values[k])
          k==1 && println("Fitting to model: ", fit_result.model_str)
          plot(xpts, fit_result.fit_curve(xpts),label="fit", color=ax[:lines][end][:get_color](), linewidth=1)
          Tvec[k] = fit_result.fit_params[fit_param_name]
          dTvec[k] = fit_result.errors[fit_param_name]
      catch
          Tvec[k] = NaN
      end
    end
  end
  label_x = collect(keys(axes))[1]
  units = data[2][group]["units"][label_x]
  if units != nothing
    label_x = string(label_x, " (", units, ")" )
  end
  xlabel(label_x)
  if cals
    ylabel(L"\langle Z\rangle")
  else
    ylabel(string(quad, "(Voltage)"))
  end
  label_y = collect(keys(axes))[2]
  units = data[2][group]["units"][label_y]
  if show_legend
    label_y = collect(keys(axes))[2]
      units = data[2][group]["units"][label_y]
    if units != nothing
      label_y = string(label_y, " (", units, ")" )
  end
    legend(title=label_y, bbox_to_anchor=[1.05,1],loc=2,borderaxespad=0)
  end
  if ~isempty(fit_name)
    plr = subplot(1,2,2)
    errorbar(ypoints_values, Tvec, dTvec, marker = "o", markersize=4)
    ylabel(fit_param_name) #TODO: unit? Proper parameter name, e.g., T1?
    xlabel(y_label)
    plr[:axis](ymin = 0)
    subplots_adjust(wspace=0.3)
  end
  title(get_partial_filename(data[2][group]["filename"]))
  ~isempty(save_fig) && savefig(string(splitext(data[2][group]["filename"])[1],'-',group,'.', save_fig), bbox_inches = "tight")
  return xpoints_values[1:length(data_values[1])], data_values, (Tvec, dTvec)
end

function plot2D_matlab(data, quad = :real; normalize=false)
  fig = figure("pyplot_surfaceplot",figsize=(5,3))
  ax = gca()
  ax[:ticklabel_format](useOffset=false)
  data_quad = eval(quad).(data["data"])
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

function pauli_set_plot(rho; rho_ideal=[], fig_width=10, fig_height=3.5, bar_width=0.6)
    pauli_vec, pauli_ops = rho2pauli(rho)
    figure(figsize=(fig_width,fig_height))
    ind = 1:length(pauli_vec)
    if ~isempty(rho_ideal)
        pauli_vec_ideal, _  = rho2pauli(rho_ideal)
        bar(ind, pauli_vec_ideal, bar_width, color="green", label=L"$\rho_{ideal}$")
    end
    bar(ind, pauli_vec[:], bar_width, label=L"$\rho_{actual}$")
    xticks(ind, map(string,pauli_ops[:]))
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

function get_partial_filename(filename, num_dirs = 3)
  cur_path = ""
  for n=1:num_dirs+1
    filename, filename_dir = splitdir(filename)
    cur_path = n==1 ? filename_dir : joinpath(filename_dir, cur_path)
  end
  return cur_path
end

function load_T1_series(datapath::AbstractString, numstart::Int, numend::Int, group, subdir=Dates.format(Dates.today(),"yymmdd");quad=:real,delta=1)
    """
      load_T1_series

    Plot multiple exponentially decaying 1D traces on top of each other from different files.
    numstart/numend: file number start/end
    group: data group name
    Optional arguments
    -------------------------
    subdir = date
    """
    T1vec = fill(NaN, numend-numstart+1)
    y0vec = fill(NaN, numend-numstart+1)
    data_temp = load_data(datapath, numstart, subdir)
    datavec = zeros(length(data_temp[1][group])-4, numend-numstart+1)
    fig = figure(figsize=(3,3))
    for (k,num) in enumerate(numstart:delta:numend)
        data = load_data(datapath, num, subdir);
        _,data_values,fit_result = Qlab.plot1D(data, group, cals=true, fit_name="t1", doplot=true, fig=fig, quad=quad)
        if fit_result != nothing
            T1vec[k] = fit_result.fit_params["T"]
            y0vec[k] = fit_result.fit_params["b"]
        end
        datavec[:,k] = data_values
    end
    return datavec, T1vec, y0vec
end

function blob(x, y, w, w_max, area; cmap=nothing)
    """
    Draws a square-shaped blob with the given area (< 1) at
    the given coordinates.
    """
    hs = sqrt(area) / 2
    xcorners = [x - hs, x + hs, x + hs, x - hs]
    ycorners = [y - hs, y - hs, y + hs, y + hs]

    plt.fill(xcorners, ycorners,
             color=cmap(round(Int64,((w + w_max) * 256 / (2 * w_max)))))
    print(((w + w_max) * 256 / (2 * w_max)))
end

function hinton(W; xlabels=nothing, ylabels=nothing, title=nothing, ax=nothing, cmap=nothing,label_top=true)

    if cmap==nothing
        cmap = ColorMap("RdBu")
    end

    if ax==nothing
        fig, ax = subplots(1, 1, figsize=(8, 6))
    else
        fig = nothing
    end

    if xlabels==nothing || ylabels==nothing
        ax.axis("off")
    end

    ax.axis("equal")
    ax.set_frame_on(false)

    height, width = size(W)

    w_max = 1.25 * maximum(abs.(Diagonal(W)))
    if w_max <= 0.0
        w_max = 1.0
    end

    ax.fill([0, width, width, 0], [0, 0, height, height],color=cmap(128))

    for x=1:width
        for y=1:height
            _x = x
            _y = y
            if real(W[x, y]) > 0.0
                blob(_x - 0.5, height - _y + 0.5, abs(W[x,
                      y]), w_max, min(1, abs(W[x, y]) / w_max), cmap=cmap)
            else
                blob(_x - 0.5, height - _y + 0.5, -abs(W[
                      x, y]), w_max, min(1, abs(W[x, y]) / w_max), cmap=cmap)
            end
        end
    end

    # color axis

    norm = matplotlib.colors.Normalize(-maximum(abs.(W)), maximum(abs.(W)))
    cax, kw = matplotlib.colorbar.make_axes(ax, shrink=0.75, pad=.1)
    matplotlib.colorbar.ColorbarBase(cax, norm=norm, cmap=cmap)

    # x axis
    ax.xaxis.set_major_locator(plt.IndexLocator(1, 0.5))

    if xlabels!=nothing
        ax.set_xticklabels(xlabels)
        if label_top!=nothing
            ax.xaxis.tick_top()
        end
    end
    ax.tick_params(axis="x", labelsize=14)

    # y axis
    ax.yaxis.set_major_locator(plt.IndexLocator(1, 0.5))
    if ylabels!=nothing
        ax.set_yticklabels((reverse(ylabels)))
    end

    ax.tick_params(axis="y", labelsize=14)

    return fig, ax
end
