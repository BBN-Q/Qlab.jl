using LsqFit

function fit_t1(xpts, ypts)
	model(t, p) = p[1]*exp(-t ./ p[2]) + p[3]
	p_guess = [maximum(ypts)-minimum(ypts), xpts[end]/3., ypts[end]]
	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	xfine = linspace(xpts[1],xpts[end],1001)
  fit_curve = (xfine, model(xfine, result.param))
  return result.param[2:3], errors[2:3], fit_curve
end

function fit_line(xpts, ypts)
	model(t, p) = p[1]*t + p[2]
	p_guess = [-ypts[indmin(abs(xpts))]/xpts[indmin(abs(ypts))], ypts[indmin(abs(xpts))]]
	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	xfine = linspace(xpts[1],xpts[end],1001)
  fit_curve = (xfine, model(xfine, result.param))
  return result.param, errors, fit_curve
end

function fit_ramsey(xpts, ypts)
	model(t, p) = p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t) + p[4]

	# Use KT estimation to get a guess for the fit
	freqs,Ts,amps = KT_estimation(ypts, xpts[2]-xpts[1], 2)
	idx = indmax(abs(amps))

	p_guess = [abs(amps[idx]), Ts[idx], freqs[idx], mean(ypts)];

	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	return result.param[2:3], errors[2:3]
end

function fit_sin(xpts, ypts)
	model(t, p) = p[1]*sin(2π*p[2] .* t) + p[3]

	# Use KT estimation to get a guess for the fit
	freqs,Ts,amps = KT_estimation(ypts, xpts[2]-xpts[1], 2)
	idx = indmax(abs(amps))

	p_guess = [amps[idx], Ts[idx], freqs[idx], mean(ypts)];

	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	return result.param, errors
end

#Fit function in McClure et al., Phys. Rev. App. 2016 (see below for details)
function fit_photon_ramsey(xpts, ypts, params)
	# input params:
	# 1 - cavity decay rate kappa (MHz)
	# 2 - detuning Delta (MHz)
	# 3 - dispersive shift 2Chi (MHz)
	# 4 - Ramsey decay time T2* (us)
	# 5 - exp(-t_meas/T1) (us), only if starting from |1>! (to include relaxation during the 1st msm't)
	# 6 - initial qubit state (0/1)
	params[1:3]*=2*pi #convert to angular frequencies
	model_0(t, p) = (-imag(exp(-(1/params[4]+params[2]*1im).*t + (p[1]-p[2]*params[3]*(1-exp(-((params[1] + params[3]*1im).*t)))/(params[1]+params[3]*1im))*1im)))
	function model(t, p)
		if params[6] == 1
			return params[5]*model_0(t, p) + (1-params[5])*model_0(t, [p[1]+pi; p[2:end]])
		else
			return model_0(t, p)
		end
	end
	p_guess = [0., 1.]
	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	xfine = linspace(xpts[1],xpts[end],1001)
	fit_curve = (xfine, model(xfine, result.param))
	return (result.param[2], errors[2], fit_curve)
end

function analyzeRB(ypts, seqlengths)
	#analyzeRB Analyzes a randomized benchmarking experiment
	#seqlengths, example: [4,8,16,32,64,128,256,384,512,640,768,896,1024,1152,1280,1408,1536,1664]

	numRepeats = length(ypts)/length(seqlengths);

	xpts = seqlengths[1 + floor((0:length(ypts)-1)./numRepeats)];

	fidelity = .5 * (1 - ypts[:]);

	model(n, p) = p[1] * (1-p[2]).^n + p[3]
	fit = curve_fit(model, xpts, fidelity, [0.5, .01, 0.5])
	xfine = linspace(seqlengths[1],seqlengths[end],1001)
	fit_curve = (xfine, model(xfine, fit.param))

	@printf("Error = %0.3f%%", fit.param[2]/2*100)
	return (xpts, fidelity, fit.param, fit_curve)
end
