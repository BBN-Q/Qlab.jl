using LsqFit

export fit_t1, fit_ramsey

function fit_t1(xpts, ypts)
	model(t, p) = p[1]*exp(-t ./ p[2]) + p[3]
	p_guess = [max(ypts)-min(ypts), xpts[end]/3., y[end]]
	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	return result.param[2:3], errors[2:3]
end

function fit_ramsey(xpts, ypts)
	model(t, p) = p[1]*exp(-t ./ p[2]).*cos(2π*p[3] .* t) + p[4]

	# Use KT estimation to get a guess for the fit
	freqs,Ts,amps = KT_estimation(ydata, xpts[2]-xpts[1], 2)
	idx = indmax(abs(amps))

	p_guess = [abs(amps[idx]), Ts[idx], freqs[idx], mean(ypts)];

	result = curve_fit(model, xpts, ypts, p_guess)
	errors = estimate_errors(result)
	return result.param[2:3], errors[2:3]
end
