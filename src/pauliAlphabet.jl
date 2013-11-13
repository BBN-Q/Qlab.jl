function pauliAlphabet(n)
	sort!(nQubitPaulis(n), by=hamming)
end

function nQubitPaulis(n)
	if n == 0
		return [""]
	end
	ops = ["I" "X" "Y" "Z"]
	paulis = ops .* nQubitPaulis(n-1)
	return paulis[:]
end

function hamming(pauliStr)
	# compute the Hamming weight of a Pauli operator (number of non-identity components)
	weight = 0
	for pauli in pauliStr
		weight += (pauli == 'I') ? 0 : 1
	end
	return weight
end