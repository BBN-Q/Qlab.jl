#Load the file loader
reload("/home/cryan/Programming/Repos/julia-tools/read_records.jl")

function filter_records(records, gCols, eCols)
	#Helper function to transform and integrate the measurement records

	#Extract the mean ground and excited traces
	meanGround = squeeze(mean(records[:,gCols,:], [2,3]), [2,3]);
	meanExcited = squeeze(mean(records[:,eCols,:], [2,3]), [2,3]);

	#Take the total mean as the shift operator
	centres = 0.5*(meanGround+meanExcited);

	#Take the difference to get the weighting and angle for the optimum quadrature
	D = meanExcited-meanGround;
	rotAngles = angle(D);

	#We square the weights to normalize to +/- 0.5
	weights = abs(D)/(sum(abs(D).^2));

	#Now we can unwind
	unwoundRecords = (records .- centres) .* exp(-1im*rotAngles);

	#Take the weighted sum and look at the good quadrature
	intRecords = squeeze(sum(weights.*real(unwoundRecords),1),1);

	return intRecords

end

function comp_to_paulis(records)
	#Helper function to convert from computational basis calibrations to pauli basis II, IZ, ZI, ZZ etc.

	#Matrix to convert from computation basis to Pauli basis diagonals
	A = 0.5*[1 1 1 1; 1 -1 1 -1; 1 1 -1 -1; 1 -1 -1 1]

	#Look at mean and variance transformed 
	meanSig = mean(records, 2)
	varSig = cov(transpose(records))

	betas = A*meanSig;
	varBetas = A*varSig*A';

	#And the standard error of the mean
	betasErrors = sqrt(diag(varBetas)/size(records,2));

	return betas, betasErrors

end

#Load the records
recordsM1 = squeeze(read_records("/home/cryan/Desktop/Records_M1"),2);
recordsM2 = squeeze(read_records("/home/cryan/Desktop/Records_M2"),2);

filteredM1 = filter_records(recordsM1, [1,2], [3,4]);
filteredM2 = filter_records(recordsM2, [1,3], [2,4]);
filteredC12 = filteredM1.*filteredM2;

betasM1, betasM1_err = comp_to_paulis(filteredM1);
betasM2, betasM2_err = comp_to_paulis(filteredM2);
betasC12, betasC12_err = comp_to_paulis(filteredC12);