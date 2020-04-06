using LsqFit
using Optim
using MicroLogging
using Roots

#Fitting a biased Lorentzian to amplitude data.
"""
Returns a five parameter, biased Lorentizian function

Args:
  x: domain of the returned function
  p: a vector of length 5 with function params --
      (amplitude, frequency, line width, linear skew, offset)
Returns:
    The computed function.

See appendix E.4 of Gao\'s thesis (2008)
"""
function bias_lorentzian(x, p)
  return p[1] ./ ((x .- p[2]).^2 .+ (p[3]/2)^2) .+ p[4]*(x .- p[2]) .+ p[5];
end

"""
Fit frequency-amplitude data to a skewed Lorentzian using nonlinear least-squares.
Model is y = a / ((x - b)^2 + (c/2)^2) + d*(x - b) + e

Args:
  freq: Frequency data.
  data: Data to be fitted.
  vars: Variance of data. Defaults to empty.
Returns:
  fit: LsqFitResult object.
  fit_function(x): Biased Lorentzian at fitted parameters.
"""
function fit_biased_lorentzian(freq, data, vars=[])

  BL_fit_dict(p) = Dict("a" => p[1], "b" => p[2], "c" => p[3],
                                "d" => p[4], "e" => p[5])
  p0 = initial_guess_blorentz(freq, data)
  return generic_fit(freq, data, bias_lorentzian, p0, BL_fit_dict,
    "y = a / ((x - b)^2 + (c/2)^2) + d*(x - b) + e", yvars=vars)
end

function initial_guess_blorentz(xpts, ypts)
  """
  Returns a length five vector with an initial guess of
  Lorentizian function parameters.
  """
  # find offset
  e = median(ypts)
  if abs(maximum(ypts) - median(ypts)) <= abs(minimum(ypts) - median(ypts))
        idx = findmin(ypts)[2]
        direction = -1
    else
        idx = argmax(ypts)
        direction = 1
  end
  # center frequency
  b = xpts[idx]
  half = abs(median(ypts) + ypts[idx]) / 2.
  if direction == 1
    idx_l = findfirst(x -> x > half, ypts)
    idx_r = findlast(x -> x > half, ypts)
  else
    idx_l = findfirst(x -> x < half, ypts)
    idx_r = findlast(x -> x < half, ypts)
  end
  c = abs(xpts[idx_l] - xpts[idx_r])
  a =  direction * c^2 * abs(ypts[idx] - e) / 4
  # slope of linear skew
  d = (ypts[end] - ypts[1])/(xpts[end] - xpts[1])
  return [a, b, c, d, e];
end

struct CircleFitParams
  """Container for the result of a circle fit.
    f0: Resonator center frequency.
    Qi: Resonator internal quality factor.
    Qc: Resonator coupling quality factor.
    ϕ: Impedance mismatch angle.
    τ: System phase delay.
    α: System gain angle.
    A: System gain amplitude.
  """
  f0
  Qi
  Qc
  ϕ
  τ
  α
  A
end

struct CircleFitResult
  fit_params::CircleFitParams
  sq_error::Float64
  SNR::Float64
  errors::CircleFitParams
  fit_curve::Function
end

#Fitting to the resonance circle of a quarter-wave resonator
function fit_resonance_circle(freq::Vector{T}, data::Vector{Complex{T}}; kwargs...) where T <: AbstractFloat
  """
  Fits complex-valued data

  Args:
    freq: Frequency data.
    data: Complex S21 data.
    kwargs: Use these if you already know something about the system and would like to override the
      following properties
      τ: Cable delay
      α: System phase shifted
      A: Overall system gain
  Returns:
    Result: The result of the fit.
  """
  @assert length(freq) == length(data) "Frequency and Data vectors must have same length!"
  @assert length(data) > 20 "Too few points."

  kwdict = Dict(kwargs)

  function χ2(data, R, xc, yc)
    x = real.(data)
    y = imag.(data)
    d = sqrt.( (x-xc).^2 + (y-yc).^2 ) - R
    return sum(d.^2)
  end

  #Fit to cable delay
  haskey(kwdict, :τ) ? τ = kwdict[:τ] : τ = fit_delay(freq, data)

  Sp = exp.(1im * 2. * π * freq * τ) .* data
  #Get best-fit circle and translate to origin.
  R, xc, yc = fit_circle(real.(Sp), imag.(Sp))

  #@assert χ2 < 10 "Could not fit circle to delay-corrected data. χ² = $(χ2)."
  @debug "Fit circle to delay corrected data with χ² = $(χ2(Sp, R, xc, yc))."

  #Translate circle to origin and fit overall phase delay and scaling
  St = Sp - (xc + 1im * yc)
  _, _, θ = fit_phase(freq, St)
  #find f → ∞ point -- this part seems fragile...
  β = mod(θ + π, π)
  P = xc + R*cos(β) + 1im*(yc + R*sin(β))

  @debug "Found scale parameters A = $(abs(P)), α = $(angle(P))."

  haskey(kwdict, :A) ? A = kwdict[:A] : A = abs(P)
  haskey(kwdict, :α) ? α = kwdict[:α] : α = angle(P)

  #Calibrate out α, A and get impedance mismatch angle
  Sc = Sp .* exp.(-1im * α) / A
  Rc, xcc, ycc = fit_circle(real.(Sc), imag.(Sc))

  #@assert χ2 < 10 "Could not fit circle to calibrated data. χ² = $(χ2)."
  @debug "Fit circle to calibrated data with χ² = $(χ2(Sc, Rc, xcc, ycc))."

  N = length(data)
  r = sqrt.((real(Sc) - xcc).^2 + (imag(Sc) - ycc).^2)
  SNR = Rc / sqrt(sum((r - Rc).^2)/(N-1))

  @debug "Data has an SNR of $(SNR)."

  ϕ = -asin(ycc/Rc)
  @debug "Found impedance mismatch angle: ϕ = $(ϕ)"
  #Final fit to phase to extract resonant frequency and Q
  Sct = Sc - (xcc + 1im * ycc)
  f0, Qr, _ = fit_phase(freq, Sct)

  Qc_cplx = Qr/(2 * Rc * exp(-1im * ϕ))
  Qi = 1/(1/Qr - real(1/Qc_cplx))
  Qc = 1/abs(1/Qc_cplx)

  @debug "Found quality factors: Qᵢ = $(Qi), Qc = $(Qc)."

  fit_params = CircleFitParams(f0, Qi, Qc, ϕ, τ, α, A)

  refined_params, errors = refine_circle_fit(freq, data, fit_params)

  χ0 = sum(abs.(data - lorentzian_resonance(fit_params, freq)).^2)
  χ1 = sum(abs.(data - lorentzian_resonance(refined_params, freq)).^2)
  @debug "Circle fit found χ² = $(χ0), refined to χ² = $(χ1)"

  if χ0 > χ1
    params = refined_params
    χ = χ1
  else
    params = fit_params
    χ = χ0
  end

  return CircleFitResult(params,
                          χ,
                          SNR,
                          errors,
                          x->lorentzian_resonance(fit_params, x))
end

function refine_circle_fit(freq, data, fit_params)
  """Refine circle fit and get errors using LsqFit.jl"""

  model(x, p) = abs.(lorentzian_resonance(p, x))

  initial_guess = [fit_params.f0, fit_params.Qi, fit_params.Qc,
                  fit_params.ϕ, fit_params.τ, fit_params.α, fit_params.A]

  fit = curve_fit(model, freq, abs.(data), initial_guess)
  errors = margin_error(fit)

  refined_params = CircleFitParams(fit.param[1], fit.param[2], fit.param[3],
                fit.param[4], fit.param[5], fit.param[6], fit.param[7])

  errs = CircleFitParams(errors[1], errors[2], errors[3],
                errors[4], errors[5], errors[6], errors[7])

  return (refined_params, errs)
end


function fit_phase(freq, data)
  """Fit phase of resonance.
  Args:
    freq: Frequency data.
    data: Complex S21 data.
  Returns:
    f0: Center frequency.
    Q: Quality factor.
    Θ0: Offset phase angle.
  """
  model(x, p) = p[1] + 2. * slope * atan.(2 * p[2] *(1. - x / p[3]))
  ϕ = unwrap(angle.(data))
  #guesses for initial parameters
  idx = findmin(abs.(ϕ - mean(ϕ)))[2]
  if mean(ϕ[1:9]) > mean(ϕ[end-9:end])
    j = findfirst(x -> x - ϕ[idx] < π/2., ϕ)
    k = findfirst(x -> x - ϕ[idx] < -π/2., ϕ)
    slope = 1
  else
    j = findfirst(x -> x - ϕ[idx] > π/2., ϕ)
    k = findfirst(x -> x - ϕ[idx] > -π/2., ϕ)
    slope = -1
  end
  Qguess = freq[idx]/abs.(freq[j] - freq[k])
  fit = curve_fit(model, freq, ϕ, [ϕ[idx], Qguess, freq[idx]])
  χ = sum(fit.resid.^2)

  @debug "Frequency vs. Phase fit found: f₀ = $(fit.param[3]), Q = $(fit.param[2]), θ₀ = $(fit.param[3])"
  @debug "Phase fit χ² = $(χ)"
  #@assert χ < 10 "Could not fit phase: χ² = $(χ)"

  return fit.param[3], fit.param[2], fit.param[1]
end

function fit_delay(freq, data)
  """Fit overall phase delay imposed by system transfer function.

  Args:
    freq: Frequency data.
    data: Complex S21 data.
  Returns:
    τ: Linear phase delay
  """
  function delay_model(x)
    ddata = data .* exp.(2. * π * 1im * freq * x)
    R, xc, yc = fit_circle(real.(ddata), imag.(ddata))
    d = sqrt.((real.(ddata) - xc).^2 + (imag.(ddata) - yc).^2) - R
    return sum(d.^2)
  end

  ϕ = unwrap(angle.(data))
  linfit(x,p) = p[1] + x * p[2]
  fit = curve_fit(linfit, freq, ϕ, [mean(ϕ), 0])
  ϕ0 = fit.param[2]
  result = optimize(delay_model, -1.5*abs(ϕ0), 1.5*abs(ϕ0)) #would prefer 1D gradient descent
  τ = Optim.minimizer(result)

  χ2 = delay_model(τ)
  @debug "Cable delay fit found: $(τ) with χ² = $(χ2)."
 #@assert χ2 < 10 "Could not calibrate out cable delay: χ² = $(χ2)."

  return τ
end


function lorentzian_resonance(p::CircleFitParams, f)
  """
  Return a resonance model in S21 amplitude over the range [x] and with
  parameters [p].

  See appendix E of Gao 2008 p. 155 Eq E.1
  """
  Qc_cplx = 1 / ((1/p.Qc)*exp(1im*p.ϕ))
  Q = 1/(1/p.Qi + real(1/Qc_cplx))
  return p.A*exp(1im*p.α).*exp.(-2π*1im*f*p.τ).*(1 - (Q/p.Qc)*exp(1im*p.ϕ)./(1 + 2*1im*Q*(f/p.f0 - 1)));
end

function lorentzian_resonance(p::Array, f)
  return lorentzian_resonance(CircleFitParams(p[1], p[2], p[3], p[4], p[5], p[6], p[7]), f)
end

function simulate_resonance(kwargs...)
  """
  Return synthetic data to test fit functions

  Returns: xpts, ypts, curve_params
  """
  f0 = 6.5 + randn()
  Qc = 120000 + randn()*2000
  Qi = 6e5 + randn()*1e4
  τ = 1.7
  ϕ = 0.17 * π
  Q = 1 ./ (1/Qi + real(1 ./ Qc*exp.(-1im*ϕ)))
  α = 0.35 * π
  A = 0.23
  df = f0/Q;
  freqs = range(f0 - 6*df, stop=f0 + 6*df, length=401)
  p = CircleFitParams(f0, Qi, Qc, ϕ, τ, α, A)
  data = lorentzian_resonance(p, freqs);
  # add some noise
  rand_phase = 2π*0.001*randn(length(freqs));
  rand_amp = 0.002*randn(length(freqs));
  data = data.*exp(-1im*rand_phase).*(1+rand_amp);

  return freqs, data, p
end

function fit_circle(x, y; refine=false)
  """Algebraic fit of (x,y) points to a circle. See:
      N. Chernov and C. J. Lesort, Least Squares Fitting of Circles,
        Journal of Mathematical Imaging and Vision, 23: 239-252, 2005.

      Args:
        x, y: X and Y points to be fitted to a circle.
      Returns:
        R: Circle radius.
        xc: x-coordinate of circle center.
        yc: y-coordinate of circle center.
  """
  @assert length(x) == length(y) "X and Y vector lengths must be equal for circle fit!"
  n = length(x)
  z = x.^2 + y.^2
  Mx = sum(x)
  My = sum(y)
  Mz = sum(z)
  Mxx = sum(x.^2)
  Myy = sum(y.^2)
  Mzz = sum(z.^2)
  Mxy = sum(x .* y)
  Mxz = sum(x .* z)
  Myz = sum(y .* z)
  M = [Mzz Mxz Myz Mz; Mxz Mxx Mxy Mx; Myz Mxy Myy My; Mz Mx My n];
  B = [0 0 0 -2; 0 1 0 0; 0 0 1 0; -2 0 0 0]

  F(η) = det(M - η*B)
  #We use Newton's method to guarantee convergence to the smallest positive eigenvalue
  ηs = newton(F, 0.0)
  ev = nullspace(M - ηs*B)

  xc = -ev[2]/2/ev[1]
  yc = -ev[3]/2/ev[1]
  R = sqrt(ev[2]^2 + ev[3]^2 - 4*ev[1]*ev[4]) / 2 ./ abs(ev[1])

  return R, xc, yc
end

function fit_circle_LM(x, y, initial_guess)
    """Iterative algorithm for refining a circle fit, using radial weighting."""

    function resid(p)
        C = sqrt.((x-p[2]).^2 + (y-p[3]).^2)
        W = 1 ./ sqrt.((p[2]-x).^2 + (p[3]-y).^2)
        return (p[1] - C) .* W
    end

    T = eltype(x)

    result = LsqFit.lmfit(resid, initial_guess, T[])

    return result.param[1], result.param[2], result.param[3]
end
