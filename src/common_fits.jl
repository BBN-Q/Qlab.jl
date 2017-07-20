# Copyright 2017 Raytheon BBN Technologies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the license at
#
#		http://www.apache.org/licenses/LICENSE-2.0

using LsqFit

## TODO: Devise a way to get specify a fit as a string...
## TODO: Fix types...

"""
FitResult type that stores the output of a fit.
"""
struct FitResult
    fit_params::Dict{String,Float64}
    sq_error::Float64
    Nσ::Float64
    errors::Dict{String,Float64}
    fit_curve::Function
    model_str::String
end

"""
`generic_fit(xpts, ypts, model, initial_guess, fit_params, model_string; yvars=[])`
Helper function to return a fit to a model.
# Arguments
* xpts: Independent fit variable.
* ypts: Dependent fit variable.
* model: Fit function(x, p) of data x and parameter vector p.
* initial_guess: Initial guess of fit parameters.
* fit_params::function(Array{Float64} => Dict{String, Float64}): Result dict of fit parameters.
* model_string::String: String indicating fit function.
* yvars: Variance of y points. Defaults to empty vector.
# Returns
* A FitResult object for the specified fit.
"""
function generic_fit(xpts, ypts, model, initial_guess, fit_params, model_string::String; yvars=[])
    @assert length(xpts) == length(ypts) "X and Y data length must match."
    if isempty(yvars)
      result = curve_fit(model, xpts, ypts, initial_guess)
      errors = estimate_errors(result)
      sq_error = sum(((model(xpts, result.param) - ypts).^2))
    else
      @assert length(ypts) == length(yvars) "Y data and Y variance lengths must match."
      result = curve_fit(model, xpts, ypts, 1./sqrt(yvars), initial_guess)
      errors = estimate_errors(result)
      sq_error = sum(((model(xpts, result.param) - ypts).^2)./yvars)
    end
    # Compute badness of fit:
    # Under the null hypothesis (that the model is valid and that the observations
    # do indeed have Gaussian statistics), the mean squared error is χ² distributed
    # with `dof` degrees of freedom. We can quantify badness-of-fit in terms of how
    # far the observed MSE is from the expected value, in units of σ = 2dof (the expected
    # standard deviation for the χ² distribution)
    dof = length(xpts)-length(initial_guess)
    Nσ = sq_error/sqrt(2dof) - dof/sqrt(2dof)
    return FitResult( fit_params(result.param),
              sq_error,
              Nσ,
              fit_params(errors),
              xpts->model(xpts,result.param),
              model_string)
end

""" `fit_t1(xpts, ypts, yvars=[])`
Fit to exponential T decay y = a exp(-t/T) + b.
"""
function fit_t1(xpts, ypts, yvars=[])
    T1_fit_dict(p) = Dict("a" => p[1], "T" => p[2], "b" => p[3])
    model(t, p) = p[1]*exp.(-t ./ p[2]) + p[3]
    p_guess = [ypts[1]-ypts[end], xpts[end]/3., ypts[end]]
    return generic_fit(xpts, ypts, model, p_guess, T1_fit_dict,
              "a * exp(-t/T) + b", yvars=yvars)
end

""" `fit_line(xpts, ypts, yvars=[])`
Fit to linear model y = ax + b.
"""
function fit_line(xpts, ypts, yvars=[])
    model(t, p) = p[1]*t + p[2]
    p_guess = [-ypts[indmin(abs.(xpts))]/xpts[indmin(abs.(ypts))], ypts[indmin(abs.(xpts))]]
    return generic_fit(xpts, ypts, model, p_guess,
    x->Dict("a" => x[1], "b"=>x[2]), "ax + b", yvars=yvars)
end

""" `fit_ramsey(xpts, ypts, yvars=[])`
Fit data to a Ramsey decay of the form a*exp(-t ./ T).*cos(2πf .* t + ϕ) + b
with optional weights for each point.
"""
function fit_ramsey(xpts, ypts, yvars=[])

    fit_dict1(p) = Dict("a"=>p[1], "T"=>p[2], "f"=>p[3], "ϕ"=>p[4], "b"=>p[5])
    model(t, p) = p[1]*exp.(-t ./ p[2]).*cos.(2π*p[3] .* t + p[4]) + p[5]
    # Use KT estimation to get a guess for the fit
    freqs,Ts,amps = KT_estimation(ypts-mean(ypts), xpts[2]-xpts[1], 1)
    p_guess = [abs(amps[1]), Ts[1], freqs[1], angle(amps[1]), mean(ypts)];

    return generic_fit(xpts, ypts, model, p_guess,
              fit_dict1, "a*exp(-t/T).*cos(2πf .* t + ϕ) + b)", yvars=yvars)
end

""" `fit_twofreq_ramsey(xpts, ypts, yvars=[])`

Fit the Ramsey data to a one-frequency model a*exp(-t ./ T).*cos(2πf .* t + ϕ) + b
or a two-frequency model
a₁*exp(-t ./ T₁).*cos(2πf₁ .* t + ϕ₁) + a₂*exp(-t ./ T₂).*cos(2πf₂ .* t + ϕ₂) + b
optionally taking variances into account, and report which should be
used based on Akaike's Information Criterion (AIC). This assumes
Gaussian statistics for each observation.
"""
function fit_twofreq_ramsey(xpts, ypts, yvars=[])
    fit_dict2(p) = Dict("a₁"=>p[1], "T₁"=>p[2], "f₁"=>p[3], "ϕ₁"=>p[4],
    "a₂"=>p[5], "T₂"=>p[6], "f₂"=>p[7], "ϕ₂"=>p[8],
    "b"=>p[9])
    model2(t, p) = ( p[1] * exp.(-t ./ p[2]) .* cos.(2π * p[3] .*t + p[4]) +
    p[5] * exp.(-t ./ p[6]) .* cos.(2π * p[7] .*t + p[8]) + p[9] )

    #Use KT estimation to get a guess for the fit
    freqs,Ts,amps = KT_estimation(ypts-mean(ypts), xpts[2]-xpts[1], 2)
    phases = angle(amps)
    amps = abs(amps)
    p_guess = [amps[1], Ts[1], freqs[1], phases[1], amps[2], Ts[2], freqs[2], phases[2], mean(ypts)]

    fitresult2 = generic_fit(xpts, ypts, model2, p_guess, fit_dict2,
    "a₁*exp(-t ./ T₁).*cos(2πf₁ .* t + ϕ₁) + a₂*exp(-t ./ T₂).*cos(2πf₂ .* t + ϕ₂) + b",
    yvars=yvars)

    fitresult1 = fit_ramsey(xpts, ypts, yvars)

    # AIC computation
    k1 = 5
    k2 = 9
    corr(k,n) = (k+1)*(k+1)/(n-k-2)
    aicc(e,k,n) = 2 * k + e + corr(k,n)

    aic = aicc(fitresult1.sq_error,k2,length(xpts)) - aicc(fitresult1.sq_error,k1,length(xpts))

    return [ fitresult1, fitresult2,
            (aic > 0) ? 1 : 2,
            aic ]
end

"""`fit_sin(xpts, ypts, yvars)`
Fit to a sine a*sin(2πf t) + b
"""
function fit_sin(xpts, ypts, yvars=[])
    model(t, p) = p[1]*sin.(2π*p[2] .* t) + p[3]
    # Use KT estimation to get a guess for the fit
    freqs,Ts,amps = KT_estimation(ypts, xpts[2]-xpts[1], 2)
    idx = indmax(abs.(amps))
    p_guess = [abs(amps[idx]), Ts[idx], freqs[idx], mean(ypts)];
    return generic_fit(xpts, ypts, model, p_guess,
              x->Dict("a"=>x[1], "f"=>x[2], "b"=>x[3]),
              "a*sin(2πf t) + b", yvars=yvars)
end

## TODO: Change this to use generic_fit
"""`fit_photon_ramsey(xpts, ypts, params)`
Fit function in McClure et al., Phys. Rev. App. 2016. Params are:
1. cavity decay rate κ (MHz)
2. detuning Δ (MHz)
3. dispersive shift 2χ (MHz)
4. Ramsey decay time T₂* (μs)
5. exp(-t_meas/T₁) (μs)
6. initial qubit state (0/1)
"""
function fit_photon_ramsey(xpts, ypts, params)
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

"""
analyzeRB(ypts, seqlengths; purity=False)
analyzeRB Analyzes a randomized benchmarking experiment
seqlengths, example: [4,8,16,32,64,128,256,384,512,640,768,896,1024,1152,1280,1408,1536,1664]
if purity = true, fit for incoherent errors. See J.J. Wallman et al., New J. Phys. 17, 113020 (2015)
Assuming data of the form <Z>s1	, <Z>s2, ...., <Z>sN, <X>s1, ..., <X>sN, ..., <Y>s1, ..., <Y>sN
with seqlengths = [s1, s2, ..., sN]
"""
function analyzeRB(ypts, seqlengths; purity=false)
    # figure out how many sequences of each length we have
    num_repeats = length(ypts) ÷ length(seqlengths) ÷ (purity ? 3 : 1)
    xpts = repeat(seqlengths, inner=num_repeats)
    if purity
        # compute length of Bloch vector
        data = vec(sum(reshape(ypts, (length(ypts) ÷ 3),3).^2,  2))
    else
        # otherwise convert <Z> to prob of 0
        data = .5 * (1 - vec(ypts))
    end
    model(n, p) = p[1] * (1-p[2]).^n + p[3]
    fit = curve_fit(model, xpts-purity, data, [0.5, .01, 0.5]) #fit to ...^(n-1) for purity
    xfine = linspace(seqlengths[1],seqlengths[end],1001)
    fit_curve = (xfine, model(xfine, fit.param))
    errors = estimate_errors(fit)
    if purity
        @printf("ϵ_inc = %0.3f%%\n", 0.5*(1-sqrt(1-fit.param[2]))*100)
    else
        @printf("ϵ = %0.3f%%\n", fit.param[2]/2*100)
    end
    return (xpts, data, fit.param, fit_curve, errors)
end
