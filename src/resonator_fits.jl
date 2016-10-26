using LsqFit

function bias_lorentz(p, x)
  """
  Returns a five parameter, biased Lorentizian function

  p : a vector of length 5 with function params --
            (amplitude, frequency, line width, linear skew, offset)
  x      : domain of the returned function

  See appendix E.4 of Gao\'s thesis (2008)
  """
  return p[1] ./ ((x - p[2]).^2 + (p[3]/2)^2) + p[4]*(x - p[2]) + p[5];
end

function lorentzian_resonance(p, f)
  """
  Return a resonance model in S21 amplitude over the range [x] and with
  parameters [p].

  See appendix E of Gao 2008 p. 155 Eq E.1
  """
  f0 = p[1]; # resonant frequency
  ϕ  = p[2]; # angle offset in the complex plane due to mismatch
  Q  = p[3]; # total quality factor
  Qc = p[4]; # coupling quality facto
  τ  = p[5]; # phase delay
  α  = p[5]; # loss
  A  = p[6]; # amplitude

  return A*exp(1im*α).*exp(-2π*1im*f*τ).*(1 - (Q/abs(Qc))*exp(1im*ϕ)./(1 + 2*1im*Q*(f/f0 - 1)));
end

function initial_guess(xpts, ypts)
  """
  Returns a length five vector with an initial guess of
  Lorentizian function parameters.
  """

  xpts = abs(xpts);
  ypts = abs(ypts);
  direction = 1;

  # find offset
  if maximum(ypts) - median(ypts) <= abs(minimum(ypts) - median(ypts))
        e = minimum(ypts);
        idx = indmax(ypts);
        direction = -1;
    else
        e = maximum(ypts);
        idx = indmin(ypts);
  end

  # center frequency
  b = xpts[idx]

  # slope of linear skew
  d = (ypts[end] - ypts[1])/(xpts[end] - xpts[1]);

  # line width
  half = abs((median(ypts) + ypts[idx])) / 2;
  if maximum(ypts) - median(ypts) <= abs(minimum(ypts) - median(ypts))
        idx_left = findfirst(x -> x > half, ypts)
        idx_right = findlast(x -> x > half, ypts)
    else
        idx_left = findfirst(x -> x < half, ypts)
        idx_right = findlast(x -> x < half, ypts)
  end
  c = abs(xpts[idx_left] - xpts[idx_right]);

  # amplitude
  a = direction * c^2 * (maximum(ypts) - minimum(ypts)) / 4;

  return [a, b, c, d, e];
end

function simulate_resonance(kwargs...)
  """
  Return synthetic data to test fit functions

  Returns: xpts, ypts, curve_params
  """

  f0 = 6.5 + randn();
  Qc = 10000 + randn()*2000;
  Qi = 5e4;
  τ = 1.74*π;
  ϕ = 2.1*π;
  Q = 1 ./ (1/Qi + real(1 ./ Qc*exp(1im*ϕ)));
  α = 1.2;
  A = 0.73;
  df = f0/Q;
  freqs = linspace(f0 - 6*df, f0 + 6*df, 401);
  data = lorentzian_resonance([f0, ϕ, Q, Qc, τ, α, A], freqs);
  # add some noise
  rand_phase = 2π*0.002*randn(length(freqs));
  rand_amp = 0.01*randn(length(freqs));
  data = data.*exp(-1im*rand_phase).*(1+rand_amp);

  return freqs, data, [f0, ϕ, Q, Qc, τ, α, A]
end
