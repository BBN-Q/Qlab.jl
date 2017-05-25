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
  e = median(ypts)
  if abs(maximum(ypts) - median(ypts)) <= abs(minimum(ypts) - median(ypts))
        idx = indmin(ypts)
        direction = -1
    else
        idx = indmax(ypts)
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
function fit_resonance_circle{T <: AbstractFloat}(freq::Vector{T}, data::Vector{Complex{T}}; kwargs...)
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

  #Fit to cable delay
  if haskey(kwdict, :τ)
    τ = kwdict[:τ]
  else
    τ = fit_delay(freq, data)
  end
  Sp = exp(1im * 2. * π * freq * τ) .* data
  #Get best-fit circle and translate to origin.
  R, xc, yc = fit_circle(real(Sp), imag(Sp))
  χ2 = sum(R^2 - (real(Sp) - xc).^2 - (imag(Sp) - yc).^2)
  @assert χ2 < 5 "Could not calibrate out cable delay: χ^2 = $χ2."
  #Translate circle to origin and fit overall phase delay and scaling
  St = Sp - (xc + 1im * yc)
  _, _, θ = fit_phase(freq, St)
  #find f → ∞ point -- this part seems fragile...
  P = xc + R*cos(θ + π) + 1im*(yc + R*sin(θ + π))
  if haskey(kwdict, :A)
    A = kwdict[:A]
  else
    A = abs(P)
  end
  if haskey(kwdict, :α)
    α = kwdict[:α]
  else
    α = angle(P)
  end
  #Calibrate out α, A and get impedance mismatch angle
  Sc = Sp .* exp(-1im * α) / A
  Rc, xcc, ycc = fit_circle(real(Sc), imag(Sc))
  ϕ = -atan2(ycc, 1 - xcc)
  #Final fit to phase to extract resonant frequency and Q
  Sct = Sc - (xcc + 1im * ycc)
  f0, Qr, _ = fit_phase(freq, Sct)
  Qc_cplx = Qr * exp(-1im * ϕ)/(2 * Rc)
  Qc = 1 / real(1 / Qc_cplx)
  Qi = 1/(1/Qr - 1/Qc)
  return CircleFitResult(f0, Qi, Qc, ϕ, τ, α, A)

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
    return freq[idx], 0, 0
  end
  Qguess = freq[idx]/abs(freq[j] - freq[k])
  fit = curve_fit(model, freq, ϕ, [ϕ[idx], Qguess, freq[idx]])
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
    ddata = data .* exp(2. * π * 1im * freq * x)
    R, xc, yc = fit_circle(real(ddata), imag(ddata))
    return sum(R.^2 - (real(ddata) - xc).^2 - (imag(ddata) - yc).^2)
  end
  ϕ = unwrap(angle(data))
  linfit(x,p) = p[1] + x * p[2]
  fit = curve_fit(linfit, freq, ϕ, [mean(ϕ), 0])
  ϕ0 = fit.param[2] / (2 * π)
  result = optimize(delay_model, -abs(ϕ0), abs(ϕ0)) #would prefer 1D gradient descent
  return Optim.minimizer(result)
end


function lorentzian_resonance(p::CircleFitResult, f)
  """
  Return a resonance model in S21 amplitude over the range [x] and with
  parameters [p].

  See appendix E of Gao 2008 p. 155 Eq E.1
  """
  Q = 1 ./ (1/p.Qi + real(1 ./ p.Qc*exp(-1im*p.ϕ)) );
  return p.A*exp(1im*p.α).*exp(-2π*1im*f*p.τ).*(1 - (Q/abs(p.Qc))*exp(1im*p.ϕ)./(1 + 2*1im*Q*(f/p.f0 - 1)));
end

function lorentzian_resonance(p::Array, f)
  return lorentzian_resonance(CircleFitResult(p[1], p[2], p[3], p[4], p[5], p[6], p[7]), f)
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
  Q = 1 ./ (1/Qi + real(1 ./ Qc*exp(-1im*ϕ)))
  α = 0.35 * π
  A = 0.23
  df = f0/Q;
  freqs = linspace(f0 - 6*df, f0 + 6*df, 401)
  p = CircleFitResult(f0, Qi, Qc, ϕ, τ, α, A)
  data = lorentzian_resonance(p, freqs);
  # add some noise
  rand_phase = 2π*0.001*randn(length(freqs));
  rand_amp = 0.002*randn(length(freqs));
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
  if abs(xc - yc) < eps()
    warn("The moment matrix of the circle fit has repeated eigenvalues; the fit has probably failed.")
  end
  return R, xc, yc
end
