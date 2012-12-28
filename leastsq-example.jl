#example use of Nonlinear fit packge for an exponential model
require("NonlinearFit")
using NonlinearFit

# produces the vector of differences between the model (p = [scale, tau]) and data
function yfit(p, xpts, data)
	return p[1]*exp(-xpts.*p[2]) - data
end

# produces the N x 2 Jacobian (note: only depends on the model, not the data)
function gfit(p, xpts, data)
	return [exp(-xpts.*p[2]) -xpts.*p[1].*exp(-xpts.*p[2])]
end

# some example data
xpts = linspace(0,10,20)
data = exp(-xpts*2) + 0.01*randn(20)
beta = [0.5, 0.5]

beta, J = leastsq(yfit, gfit, beta, (xpts, data))

println("Found beta: $beta")