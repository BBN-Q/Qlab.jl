module NonlinearFit

require("Options")
using Base
using OptionsMod

export leastsq

function curve_fit(f::Function, xdata, ydata, p0...)
	# assumes f(xdata, params...) = ydata + epsilon
	# minimizes sum(ydata - f(xdata)).^2 using leastsq()
end

leastsq(f::Function, g::Function, x0, fargs) = leastsq(f::Function, g::Function, x0, fargs, @options)

function leastsq(f::Function, g::Function, x0, fargs, opts::Options)
	# finds min_x sum(f(x).^2) using the Levenberg-Marquardt algorithm
	# The function f should take an input vector of length n and return an output vector of length m
	# The function g is the Jacobian of f, and should be an m x n matrix
	# x0 is an initial guess for the solution
	# fargs is a tuple of additional arguments to pass to f
	# available options:
	#   tolX - search tolerance in x
	#   tolG - search tolerance in gradient
	#   maxIter - maximum number of iterations
	#   lambda - (inverse of) initial trust region radius
	n = size(x0,1)
	@defaults opts tolX=1e-8 tolG=1e-6 maxIter=100*n lambda=100.0

	# other constants
	const MAX_LAMBDA = 1e16 # maximum trust region radius
	const MIN_LAMBDA = 1e-32 # minimum trust region radius
	const MIN_STEP_QUALITY = 1e-3
	# const MIN_DIAGONAL = 1e6 # lower bound on values of diagonal matrix used to regularize the trust region step
	# const MAX_DIAGONAL = 1e32 # lower bound on values of diagonal matrix used to regularize the trust region step

	iterCt = 1
	x = x0
	delta_x = copy(x0)

	fcur = f(x, fargs...)
	residual = sse(fcur)
	println("Initial residual: $residual")

	while ( iterCt < maxIter && norm(delta_x) > tolX*(tolX + norm(x)) )
		println("x: $x")
		J = g(x, fargs...)
		# we want to solve:
		#    argmin 0.5*||J(x)*delta_x + f(x)||^2 + lambda*||J^T*J*delta_x||^2
		# We can combine the two terms by concatenating J^T*J onto the bottom of J,
		# and padding f(x) with zeros (see Ceres Solver user guide, section 13, p. 55)
		# Also use the equivalence: diagm(J.'*J) = diagm(sum(J.^2, 1))
		# Solving for the minimum gives:
		#    J^T * J_aug * delta_x == -J^T * f(x), where J_aug = [J, diagm(sqrt(lambda)*sum(J.^2,1))]
		# Again, referring to the Ceres guide, p. 58, it looks like we can drop the J^T from both sides
		# in the resulting least-squares QR problem. So, we have:
		delta_x = [J, diagm(sqrt(lambda)*sum(J.^2,1))] \ [-fcur, zeros(n)]
		println("delta_x: $delta_x")
		# if the linear assumption is valid, our new residual should be:
		predicted_residual = sse(J*delta_x + fcur)
		# check for numerical problems in solving for delta_x by ensuring that the predicted residual is smaller
		# than the current residual
		if predicted_residual > residual
			error("Error solving for delta_x: predicted residual increase.")
		end
		# try the step and compute its quality
		trial_f = f(x + delta_x, fargs...)
		trial_residual = sse(trial_f)
		println("Trial residual: $trial_residual")
		# step quality = residual change / predicted residual change
		rho = (trial_residual - residual) / (predicted_residual - residual)
		println("Step quality rho: $rho")

		if rho > MIN_STEP_QUALITY
			x += delta_x
			fcur = trial_f
			residual = trial_residual
			# increase trust region radius
			lambda = max(0.1*lambda, MIN_LAMBDA)
		else
			# decrease trust region radius
			lambda = min(10*lambda, MAX_LAMBDA)
		end
		iterCt += 1

		# should check for the following stopping conditions:
		# 1. Small gradient: norm(J^T * fcur) < tolG
		# 2. Small step size: norm(delta_x) < tolX
		# 3. iterCt > maxIter
		if norm(J.' * fcur) < tolG
			println("Stopping: small gradient.")
			break
		end
	end

	# give the user info about the stopping condition
	if iterCt >= maxIter
		println("Exceeded maximum number of iterations")
	elseif norm(delta_x) <= tolX * (tolX + norm(x))
		println("Step size too small.")
	end

	J = g(x, fargs...)

	return x, J

end

sse(x) = (x.'*x)[1]

end