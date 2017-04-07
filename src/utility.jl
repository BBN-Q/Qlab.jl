function unwrap!{T <: AbstractFloat}(Θ::Array{T}; discont=π)
  """Unwrap angles by changing deltas between values to 2π complement.

  Unwrap radian phase Θ by changing absolute jumps greater than `discont` to
  their 2π complement.

  Args:
    Θ: Input array.
    discont: Maximum discontinuity between values. Optional, default is π.

  Returns:
    out: Output array of unwrapped phases.

  See also:
    numpy.unwrap()
  """
  Δ = diff(Θ)
  Δmod = mod(Δ + π, 2 * π) - π
  Δmod[(Δmod .== π) & (Δ .> 0)] = π
  Θcorr = Δmod - Δ
  Θcorr[abs(Δ) .< discont] = 0
  return Θ .+ vcat(0, cumsum(Θcorr))
end

function unwrap{T <: AbstractFloat}(Θ::Array{T}, args...; kwargs...)
  """In-place version of unwrap."""
  unwrap!(copy(Θ), args...; kwargs...)
end
