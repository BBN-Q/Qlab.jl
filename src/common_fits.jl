using LsqFit, Distributions

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

"""
	fit_ramsey(xpts, ypts, weights=[])

Fit data to a Ramsey decay of the form p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t + p[4]) + p[5],
with optional weights for each point.
"""
function fit_ramsey(xpts, ypts)
	model(t, p) = p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t) + p[4]

	# Use KT estimation to get a guess for the fit
	freqs,Ts,amps = KT_estimation(ypts, xpts[2]-xpts[1], 2)
	idx = indmax(abs(amps))

	p_guess = [abs(amps[idx]), Ts[idx], freqs[idx], mean(ypts)];

	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	xfine = linspace(xpts[1],xpts[end],1001)
	fit_curve = (xfine, model(xfine, result.param))
	return result.param[2:3], errors[2:3], fit_curve
end

"""
	fit_twofreq_ramsey(xpts, ypts, weights, compare=false, verbose=false)

Fit a two-frequency decaying cosine to Ramsey data, with optional weights.
If `compare`, will compare the fit to a one-frequency Ramsey fit using
a statistical F-test. If `verbose`, will print results of fits.
"""
function fit_twofreq_ramsey(xpts, ypts, weights=[])
    model(t, p) = ( p[1] * exp(-t ./ p[2]) .* cos(2π * p[3] .*t + p[4]) +
		    p[5] * exp(-t ./ p[6]) .* cos(2π * p[7] .*t + p[8]) + p[9] )
    #Use KT estimation to get a guess for the fit
    freqs,Ts,amps = KT_estimation(ypts, xpts[2]-xpts[1], 2)
    inds = find(x->(x > 0), Ts)
    if length(inds) < 2
	@printf("KT Estimation for two frequency failed.\n")
	return [], []
    end
    freqs = freqs[inds]
    Ts = Ts[inds]
    amps = amps[inds]
    phases = angle(amps)
    amps = abs(amps)
    p_guess = [amps[1], Ts[1], freqs[1], phases[1], amps[2], Ts[2], freqs[2],phases[2], mean(ypts)]
    if isempty(weights)
  	result2 = curve_fit(model, xpts, ypts, p_guess)
    else
	result2 = curve_fit(model, xpts, ypts, weights, p_guess)
    end
    errors2 = estimate_errors(result2)
    # compute weighted least squares error, assuming weights = pointwise inverse variance
    sq_err2 = ((model(xpts, result2.param) - ypts).^2).*weights

    xfine = linspace(xpts[1],xpts[end],1001)
    fit_curve2 = (xfine, model(xfine, result2.param))
    
    model1(t,p) = p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t + p[4]) + p[5]
    freqs,Ts,amps = KT_estimation(ypts, xpts[2]-xpts[1], 2)
    idx = indmax(abs(amps))
    p_guess = [abs(amps[idx]), Ts[idx], freqs[idx], angle(amps[idx]), mean(ypts)]
    if isempty(weights)
	result1 = curve_fit(model1, xpts, ypts, p_guess)
    else
	result1 = curve_fit(model1, xpts, ypts, weights, p_guess)
    end
    errors1 = estimate_errors(result1)
    # compute weighted least squares error, assuming weights = pointwise inverse variance
    sq_err1 = ((model1(xpts, result1.param) - ypts).^2).*weights

    fit_curve1 = (xfine, model1(xfine, result1.param))

    # AIC computation
    k1 = 3
    k2 = 6 
    corr(k,n) = (k+1)*(k+1)/(n-k-2)
    aicc(e,k,n) = 2 * k + e + corr(k,n)

    aic = aicc(sq_err2,k2,length(xpts)) - aicc(sq_err1,k1,length(xpts))

    return (result1.param, errors1, fit_curve1), (result2.param, errors2, fit_curve2), aic > 0 ? 1 : 2
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
