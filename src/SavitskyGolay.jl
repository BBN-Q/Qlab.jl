#Polynomial smoothing with the Savitsky Golay filters
#
# Sources
# ---------
# Theory: http://www.ece.rutgers.edu/~orfanidi/intro2sp/orfanidis-i2sp.pdf
# Python Example: http://wiki.scipy.org/Cookbook/SavitzkyGolay

function savitsky_golay(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)

#Some error checking
@assert isodd(windowSize) "Window size must be an odd integer."
@assert polyOrder < windowSize "Polynomial order must me less than window size."

halfWindow = int((windowSize-1)/2)

#Setup the S matrix of basis vectors. 
S = zeros(windowSize, polyOrder+1)
for ct = 0:polyOrder
	S[:,ct+1] = [-halfWindow:halfWindow].^(ct)
end

#Compute the filter coefficients for all orders
#From the scipy code it seems pinv(S) and taking rows should be enough
G = S*pinv(S'*S)

#Slice out the derivative order we want
filterCoeffs = G[:,deriv+1] * factorial(deriv);

#Pad the signal with the endpoints and convolve with filter 
paddedX = [x[1]*ones(halfWindow), x, x[end]*ones(halfWindow)]
y = conv(filterCoeffs, paddedX)

#Return the valid midsection
return y[2*halfWindow+1:end-2*halfWindow]

end
