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

immutable FitResult
    fit_params::Vector{Float64}
    sq_errror::Float64
    Nσ::Float64
    errors::Vector{Float64}
    fit_curve::Function
end

"""
        fit_ramsey(xpts, ypts, yvars=[])

Fit data to a Ramsey decay of the form p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t + p[4]) + p[5],
with optional weights for each point.
"""
function fit_ramsey(xpts, ypts, yvars=[])
    model(t, p) = p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t + p[4]) + p[5]

    # Use KT estimation to get a guess for the fit
    freqs,Ts,amps = KT_estimation(ypts-mean(ypts), xpts[2]-xpts[1], 1)

    p_guess = [abs(amps[1]), Ts[1], freqs[1], angle(amps[1]), mean(ypts)];

    if isempty(yvars)
        result = curve_fit(model, xpts, ypts, p_guess)
    else
        result = curve_fit(model, xpts, ypts, 1./sqrt(yvars), p_guess)
    end
    errors = estimate_errors(result)
    sq_error = sum(((model(xpts, result.param) - ypts).^2)./yvars)

    dof = length(xpts)-5
    Nσ = sq_error/sqrt(2dof) - dof/sqrt(2dof)

    xfine = linspace(xpts[1],xpts[end],1001)
    fit_curve = (xfine, model(xfine, result.param))

    return FitResult( result.param, sq_error, Nσ, errors, x->model(x,result.param) )
end

"""
        fit_twofreq_ramsey(xpts, ypts, yvars=[])

Fit one- and two-frequency decaying cosine to Ramsey data, optionally
taking variances into account, and report which should be used based
on Akaike's Information Criterion (AIC). This assumes Gaussian
statistics for each observation.
"""
function fit_twofreq_ramsey(xpts, ypts, yvars=[])
    model2(t, p) = ( p[1] * exp(-t ./ p[2]) .* cos(2π * p[3] .*t + p[4]) +
                    p[5] * exp(-t ./ p[6]) .* cos(2π * p[7] .*t + p[8]) + p[9] )
    #Use KT estimation to get a guess for the fit
    freqs,Ts,amps = KT_estimation(ypts-mean(ypts), xpts[2]-xpts[1], 2)

    phases = angle(amps)
    amps = abs(amps)
    p_guess = [amps[1], Ts[1], freqs[1], phases[1], amps[2], Ts[2], freqs[2], phases[2], mean(ypts)]

    if isempty(yvars)
        result2 = curve_fit(model2, xpts, ypts, p_guess)
        sq_err2 = sum((model2(xpts, result2.param) - ypts).^2)
    else
        result2 = curve_fit(model2, xpts, ypts, 1./sqrt(yvars), p_guess)
        sq_err2 = sum(((model2(xpts, result2.param) - ypts).^2)./yvars)
    end
    errors2 = estimate_errors(result2)
    # Compute badness of fit:
    # Under the null hypothesis (that the model is valid and that the observations
    # do indeed have Gaussian statistics), the mean squared error is χ² distributed
    # with `dof2` degrees of freedom. We can quantify badness-of-fit in terms of how
    # far the observed MSE is from the expected value, in units of σ = 2dof (the expected
    # standard deviation for the χ² distribution)
    dof2 = length(xpts)-9
    Nσ2 = (sq_err2 - dof2)/sqrt(2dof2)

    # Now do 1 freq fit for comparisson TODO: just call fit Ramsey
    model1(t,p) = p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t + p[4]) + p[5]
    freqs,Ts,amps = KT_estimation(ypts-mean(ypts), xpts[2]-xpts[1], 1)

    p_guess = [abs(amps[1]), Ts[1], freqs[1], angle(amps[1]), mean(ypts)]
    if isempty(yvars)
        result1 = curve_fit(model1, xpts, ypts, p_guess)
        sq_err1 = sum((model1(xpts, result1.param) - ypts).^2)
    else
        result1 = curve_fit(model1, xpts, ypts, 1./sqrt(yvars), p_guess)
        sq_err1 = sum(((model1(xpts, result1.param) - ypts).^2)./yvars)
    end
    errors1 = estimate_errors(result1)
    # Compute badness of fit for this model
    dof1 = length(xpts)-5
    Nσ1 = (sq_err1 - dof1)/sqrt(2dof1)

    # AIC computation
    k1 = 5
    k2 = 9 
    corr(k,n) = (k+1)*(k+1)/(n-k-2)
    aicc(e,k,n) = 2 * k + e + corr(k,n)

    aic = aicc(sq_err2,k2,length(xpts)) - aicc(sq_err1,k1,length(xpts))

    return [ FitResult( result1.param, sq_error1, Nσ1, errors1, x->model1(x,result1.param) ),
             FitResult( result2.param, sq_error2, Nσ2, errors2, x->model2(x,result2.param) )],
           (aic > 0) ? 1 : 2, 
           aic
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
