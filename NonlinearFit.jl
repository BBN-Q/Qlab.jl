Module NonlinearFit

require("Options")
using Base
using OptionsMod

function curve_fit(f, xdata, ydata, p0...)
    # assumes f(xdata, params...) = ydata + epsilon
    # minimizes sum(ydata - f(xdata)).^2 using leastsq()
end

function leastsq(f, x0, fargs, opts::Options)
    # finds min_x sum(f(x).^2) using the Levenberg-Marquardt algorithm
    # The function f may return one value or two. If f returns two values, the
    # second is interpreted as the Jacobian
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
    prevNorm = 0.0

    while (iterCt < maxIter && funEvals < maxFunEvals && delta_x > tolX)
        (fcur, J) = f(x, fargs...)
        fnorm = norm(fcur, 1)
        funEvals += 1
        if funEvals > maxFunEvals
            print("Exceeded maximum number of function evaluations, exiting...")
            break
        end
        if iterCt > 1 && abs(fnorm - prevNorm) < tolF
            print("Function decreased by less than tolF, exiting...")
            break
        end
        if fnorm < prevNorm
            # try to reduce lambda
            lambda *= max(.1, eps())
        else
            # try to increase lambda to get going in the right direction again
            while (prevNorm < from)
                lambda *= 10
                delta_x = (J.'*J + lambda*diagm(J.'*J)) \ (J.' * fnorm)
                x += delta_x
                (fcur, J) = f(x, fargs...)
                fnorm = norm(fcur, 1)
                funEvals += 1
                if funEvals > maxFunEvals
                    print("Exceeded maximum number of function evaluations, exiting...")
                    break
                end
            end
        end
        
        # find delta_x by solving the system:
        # (J^T*J + lambda * diag(J^T*J)) * delta_x = J^T * fnorm
        # also try diag(J^T*J) = sum(J.^2,1)
        delta_x = (J.'*J + lambda*diag(J.'*J)) \ (J.' * fnorm)
        
        x += delta_x
        prevNorm = fnorm
        iterCt += 1
    end
end

end