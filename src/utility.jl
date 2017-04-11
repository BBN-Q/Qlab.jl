function unwrap!{T <: AbstractFloat}(ϕ::Array{T}; discont=π)
  """Unwrap angles by changing deltas between values to 2π complement.

  Unwrap radian phase ϕ by changing absolute jumps greater than `discont` to
  their 2π complement.

  Args:
    ϕ: Input array.
    discont: Maximum discontinuity between values. Optional, default is π.

  Returns:
    out: Output array of unwrapped phases.

  See also:
    numpy.unwrap()
  """
  Δ = diff(ϕ)
  Δmod = mod(Δ + π, 2 * π) - π
  Δmod[(Δmod .== -π) & (Δ .> 0)] = π
  ϕcorr = Δmod - Δ
  ϕcorr[abs(Δ) .< discont] = 0
  return ϕ .+ vcat(0, cumsum(ϕcorr))
end

function unwrap{T <: AbstractFloat}(ϕ::Array{T}; discont=π)
  """In-place version of unwrap."""
  return unwrap!(copy(ϕ), discont=discont)
end
