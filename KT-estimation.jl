function hilbert(signal)
	# construct the Hilbert transform of the signal via the FFT
	# in essense, we just want to set negative frequency components to zero
	spectrum = fft(signal)
	spectrum[end/2+1:end] = 0
	ifft(2*spectrum)
end

function KT_estimation(data, timeStep, order)
	#Perform a modified KT estimation to obtain the estimates of the parameters 
	#from a set of data. 
	#
	# function [freqs, Tcs, amps] = KT_estimation(data, timeStep, order)
	# 
	# See ? Van Huffel, S. (1993). Enhanced resolution based on minimum variance estimation and exponential data modeling.
	#       Signal Processing, 33(3), 333-355. doi:10.1016/0165-1684(93)90130-3

	analyticSig = hilbert(data)

	#Create the raw Hankel matrix
	N = length(analyticSig)
	K = order
	M = itrunc(N/2)-1
	L = N-M+1
	H = zeros(L,M)
	for ct = 1:M
	    H[:,ct] = analyticSig[ct:ct+L-1]
	end

	#Try and seperate the signal and noise subspace via the svd
	U,S,V = svd(H)

	#Reconstruct the approximate Hankel matrix with the first K singular values
	#Here we can iterate and modify the singular values
	S_k = diagm(S[1:K])
	#Estimate the variance from the rest of the singular values
	varEst = (1/((M-K)*L)) * sum(S[K+1:end].^2)
	Sfilt = (S_k.^2 - L*varEst*eye(K)) / S_k
	Hbar = U[:,1:K] * Sfilt * V[:,1:K]'

	#Reconstruct the data from the averaged anti-diagonals
	cleanedData = zeros(N)
	tmpMat = fliplr(Hbar)
	idx = -L+1
	for ct = N:-1:1
	    cleanedData[ct] = mean(diag(tmpMat,idx))
	    idx = idx+1
	end

	#Create a cleaned Hankel matrix
	#(BRJ note: this is the same as lines 26-28. Do we want to construct a cleaned analytic signal first?)
	cleanedH = zeros(L,M)
	for ct = 1:M
	    cleanedH[:,ct] = analyticSig[ct:ct+L-1]
	end

	#Compute Q with total least squares
	#U_K1*Q = U_K2
	U,_,_ = svd(cleanedH)
	U_K = U[:,1:K]
	tmpMat = [U_K[1:end-1,:] U_K[2:end,:]]
	_,_,V = svd(tmpMat)
	n = size(U_K,2)
	V_AB = V[1:n,1+n:end]
	V_BB = V[1+n:end,1+n:end]
	Q = -V_AB/V_BB

	#Now poles are eigenvalues of Q
	poles = eig(Q, false)

	#Take the log and return the decay constant and frequency
	freqs = zeros(K)
	Tcs = zeros(K)
	for ct = 1:K
	    sk = log(poles[ct])
	    freqs[ct] = imag(sk)/2/pi/timeStep
	    Tcs[ct] = -1/real(sk)*timeStep
	end

	#Refit the data to get the amplitude
	A = zeros(N,K)
	for ct = 1:K
	    A[:,ct] = poles[ct].^(0:N-1)
	end

	amps = A\cleanedData

	return freqs, Tcs, amps

end
