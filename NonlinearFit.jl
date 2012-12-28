Module NonlinearFit

require("Options")
using Base
using OptionsMod

function curve_fit(f::Function, xdata, ydata, p0...)
    # assumes f(xdata, params...) = ydata + epsilon
    # minimizes sum(ydata - f(xdata)).^2 using leastsq()
end

function leastsq(f::Function, g::Function, x0, fargs, opts::Options)
    # finds min_x sum(f(x).^2) using the Levenberg-Marquardt algorithm
    # The function f should take an input vector of length n and return an output vector of length m
    # The function g is the Jacobian of f, and should be an m x n matrix
    # x0 is an initial guess for the solution
    # fargs is a tuple of additional arguments to pass to f
    # available options:
    #   tolX - search tolerance in x
    #   tolF - search tolerance in the function f
    #   maxFunEvals - maximum number of function evaluations
    #   maxIter - maximum number of iterations
    #   lambda - controls 
    n = len(x0)
    @defaults opts tolX=2*eps() tolF=2*eps() maxFunEvals=200*(n+1) maxIter=100*(n+1) lambda=100.0
    
    iterCt = 1
    x = x0
    delta_x = Inf
    prevNorm = Inf

    fcur = f(x, fargs...)
    residual = norm(fcur)
    funEvals = 1

    while (iterCt < maxIter && funEvals < maxFunEvals && delta_x > tolX)
        prevNorm = residual
        fcur = f(x, fargs...)
        J = g(x, fargs...)
        residual = norm(fcur)
        funEvals += 1
        if funEvals > maxFunEvals
            error("Exceeded maximum number of function evaluations, exiting...")
        end
        if iterCt > 1 && abs(residual - prevNorm) < tolF
            error("Function decreased by less than tolF, exiting...")
        end
        if residual < prevNorm
            # try to reduce lambda
            lambda *= max(.1, eps())
        else
            # try to increase lambda to get going in the right direction again
            while (residual > prevNorm)
                lambda *= 10
                if lambda > MAX_LAMBDA
                    error("Exceeded maximum trust region radius")
                end
                # use the equivalence: diagm(J.'*J) = diagm(sum(J.^2, 1))
                # delta_x = (J.'*J + diagm(lambda*sum(J.^2,1))) \ -(J.' * residual)
                # we use an additional trick... since we only need the norm to compute 
                # the above line, we can combine the vectors. i.e.:
                delta_x = [J, diagm(sqrt(lambda*sum(J.^2,1)))] \ [residual, zeros(size(x))]
                # (see Ceres section 13, p. 55 and p. 58)
                trial_x += delta_x
                fcur = f(trial_x, fargs...)
                residual = norm(fcur)
                funEvals += 1
                if funEvals > maxFunEvals
                    error("Exceeded maximum number of function evaluations, exiting...")
                end
            end
        end
        
        # find delta_x by solving the system:
        # (J^T*J + lambda * diagm(J^T*J)) * delta_x = -J^T * residual
        # see notes above
        delta_x = [J, diagm(sqrt(lambda*sum(J.^2,1)))] \ [residual, zeros(size(x))]
        
        x += delta_x
        prevNorm = residual
        iterCt += 1
    end
end

end