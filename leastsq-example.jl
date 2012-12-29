#example use of Nonlinear fit packge for an exponential model
require("NonlinearFit")
using NonlinearFit

model(xpts, p) = p[1]*exp(-xpts.*p[2])

# produces the vector of differences between the model (p = [scale, tau]) and data
f(p, xpts) = model(xpts, p) - data

# produces the N x 2 Jacobian (note: only depends on the model, not the data)
g(p, xpts) = [exp(-xpts.*p[2]) -xpts.*p[1].*exp(-xpts.*p[2])]

# some example data
xpts = linspace(0,10,20)
data = exp(-xpts*2) + 0.01*randn(20)

beta, J = leastsq(f, g, [0.5, 0.5], (xpts,))
println("Called leastsq with analytic Jacobian")
println("Found beta: $beta")

# or we can let curve_fit construct these functions for us
beta, r, J = curve_fit(model, xpts, data, [0.5, 0.5])

println("Called curve_fit")
println("Found beta: $beta")

# can also get error estimates on the fit parameters
errors = estimate_errors(beta, r, J)

println("Estimated errors: $errors")