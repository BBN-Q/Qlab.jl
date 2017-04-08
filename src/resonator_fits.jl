using LsqFit
using Optim

#Fitting a biased Lorentzian to amplitude data.
function bias_lorentzian(x, p)
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
  return p[1] ./ ((x - p[2]).^2 + (p[3]/2)^2) + p[4]*(x - p[2]) + p[5];
end

function fit_biased_lorentzian(freq, data)
  """
  Fit frequency-amplitude data to a skewed Lorentzian using nonlinear least-squares.

  Args:
    freq: Frequency data.
    data: Data to be fitted.
    alpha: Confidence limit for error calculation. Default 0.95.
  Returns:
    fit: LsqFitResult object.
    fit_function(x): Biased Lorentzian at fitted parameters.
  """
  p0 = initial_guess(freq, data)
  fit = curve_fit(bias_lorentzian, freq, data, p0)
  errors = estimate_errors(fit, 0.95)
  fit_function(x) = bias_lorentzian(x, fit.param)
  return fit, fit_function
end

function initial_guess(xpts, ypts)
  """
  Returns a length five vector with an initial guess of
  Lorentizian function parameters.
  """
  # find offset
  if abs(maximum(ypts) - median(ypts)) <= abs(minimum(ypts) - median(ypts))
        e = median(ypts)
        idx = indmin(ypts)
        direction = -1
    else
        e = median(ypts)
        idx = indmax(ypts)
        direction = 1
  end
  # center frequency
  b = xpts[idx]
  half = direction * abs(median(ypts) - ypts[idx]) / 2.
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

immutable CircleFitResult
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

#Fitting to the resonance circle of a quarter-wave resonator
function fit_resonance_circle{T <: AbstractFloat}(freq::Array{T, 1}, data::Array{Complex{T}, 1})
  """
  Fits complex-valued data

  Args:
    freq: Frequency data.
    data: Complex S21 data.
  Returns:
    Result: The result of the fit.
  """
  @assert length(freq) == length(data) "Frequency and Data vectors must have same length!"
  @assert length(data) > 20 "Too few points."

  τ, α, A = calibrate_resonance_circle(freq, data)
  scaled_data = apply_calibration(τ, α, A, freq, data)

  R, xc, yc = fit_circle(real(scaled_data), imag(scaled_data))
  ϕ = -asin(yc / R)
  td = (real(scaled_data) - xc) + 1im*(imag(scaled_data) - yc)
  f0, Qfit, _ = fit_phase(freq, td)
  Qc = Qfit / (2. * R * exp(-1im * ϕ))
  Qi = 1./ (1./Qfit - real(1./Qc))
  Qc = abs(Qc)

  return CircleFitResult(f0, Qi, Qc, ϕ, -τ, -α, A)
end

function calibrate_resonance_circle(freq, data)
  """Calibrates out cable delay and overall system gain to translate resonance
  circle to "canonical" position for fit.

  Args:
    freq: Frequency data.
    data: Complex S21 data.
  Returns:
    τ: Fitted cable delay.
    α: Overall system gain angle.
    A: Overall system gain amplitude.
  """
  #Fit to cable delay.
  τ = fit_delay(freq, data)
  Sp = exp(-1im * 2. * π * freq * τ) .* data
  #Get best-fit circle and translate to origin.
  R, xc, yc = fit_circle(real(Sp), imag(Sp))
  Strans = (real(Sp) - xc) + 1im * (imag(Sp) - yc)
  #Calculate invariant ω → ∞ point
  _, _, θ = fit_phase(freq, Strans)
  P = xc + R * cos(θ + π) + 1im * (yc + R * sin(θ + π))
  A = abs(P)
  α = angle(P)
  return τ, α, A
end

function apply_calibration(τ, α, A, freq, data)
  """Apply the calibration to the resonance data to move it to the
  canonical position.

  Args:
    τ: Fitted cable delay.
    α: Overall system gain angle.
    A: Overall system gain amplitude.
    freq: Frequency data.
    data: Complex S21 data.
  Returns:
    scaled_data: Data shifted to the canonical position for a circle fit.
  """
  data = data .* exp(-1im * 2 * π * τ * freq)
  rot = [cos(-α) -sin(-α); sin(-α) cos(-α)]
  scaled_data = zeros(Complex128, size(data))
  for j = 1:length(data)
    v = rot * [real(data[j]); imag(data[j])]
    scaled_data[j] = (v[1] + 1im * v[2]) / A
  end
  return scaled_data
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
  model(x, p) = p[1] + 2. * slope * atan(2 * p[2] *(1. - x / p[3]))
  ϕ = unwrap(angle(data))
  #first some initial guesses
  idx = indmin(abs(ϕ - mean(ϕ)))
  if mean(ϕ[1:9]) > mean(ϕ[end-9:end])
    j = findfirst(x -> x - ϕ[idx] < π/2., ϕ)
    k = findfirst(x -> x - ϕ[idx] < -π/2., ϕ)
    slope = 1
  else
    j = findfirst(x -> x - ϕ[idx] > π/2., ϕ)
    k = findfirst(x -> x - ϕ[idx] > -π/2., ϕ)
    slope = -1
  end
  #If the initial complex S21 data is very close to a circle, there is no need to
  #calibrate out cable delay so this function will fail. Instead just return a "naive"
  #guess and move on.
  if j == 0 || k == 0
    fitfunc(x) = model(x, [0, 0, freq[idx]])
    return freq[idx], 0, 0, fitfunc
  end
  Qguess = freq[idx]/abs(freq[j] - freq[k])
  fit = curve_fit(model, freq, ϕ, [ϕ[idx], Qguess, freq[idx]])
  fitfunc(x) = model(x, fit.param)
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
    data = data .* exp(-2. * π * 1im * freq * x)
    R, xc, yc = fit_circle(real(data), imag(data))
    return sum(R.^2 - (real(data) - xc).^2 - (imag(data) - yc).^2)
  end
  ϕ = unwrap(angle(data))
  linfit(x,p) = p[1] + x * p[2]
  fit = curve_fit(linfit, cat(1, freq[1:9], freq[end-9:end]), cat(1, ϕ[1:9], ϕ[end-9:end]), [mean(ϕ), 0])
  result = optimize(delay_model, -abs(fit.param[2]), abs(fit.param[2])) #would prefer 1D gradient descent
  return Optim.minimizer(result)
end


function lorentzian_resonance(p::CircleFitResult, f)
  """
  Return a resonance model in S21 amplitude over the range [x] and with
  parameters [p].

  See appendix E of Gao 2008 p. 155 Eq E.1
  """
  Q = 1 ./ (1/p.Qi + real(1 ./ p.Qc*exp(1im*p.ϕ)));
  return p.A*exp(1im*p.α).*exp(-2π*1im*f*p.τ).*(1 - (Q/abs(p.Qc))*exp(1im*p.ϕ)./(1 + 2*1im*Q*(f/p.f0 - 1)));
end

function lorentzian_resonance(p::Array, f)
  return lorentzian_resonance(CircleFitResult(p[1], p[2], p[3], p[4], p[5], p[6]), f)
end

function simulate_resonance(kwargs...)
  """
  Return synthetic data to test fit functions

  Returns: xpts, ypts, curve_params
  """
  f0 = 6.5 + randn();
  Qc = 26000 + randn()*2000;
  Qi = 1.3e5;
  τ = 1.74*π;
  ϕ = 2.1*π;
  Q = 1 ./ (1/Qi + real(1 ./ Qc*exp(1im*ϕ)));
  α = 1.2;
  A = 0.73;
  df = f0/Q;
  freqs = linspace(f0 - 6*df, f0 + 6*df, 401);
  p = CircleFitResult(f0, Qi, Qc, ϕ, τ, α, A)
  data = lorentzian_resonance(p, freqs);
  # add some noise
  rand_phase = 2π*0.002*randn(length(freqs));
  rand_amp = 0.01*randn(length(freqs));
  data = data.*exp(-1im*rand_phase).*(1+rand_amp);

  return freqs, data, p
end

function fit_circle(x, y)
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
  D, V = eig(M, B)
  D[D .< eps()] = NaN
  ev = V[:, indmin(D)]
  xc = -ev[2]/2/ev[1]
  yc = -ev[3]/2/ev[1]
  R = sqrt(ev[2]^2 + ev[3]^2 - 4*ev[1]*ev[4])/2./abs(ev[1])
  return R, xc, yc
end
