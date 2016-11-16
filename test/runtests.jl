using Qlab

function test_comp_to_paulis()
	#Load the records
	recordsM1 = squeeze(read_records("/home/cryan/Desktop/Records_M1"),2);
	recordsM2 = squeeze(read_records("/home/cryan/Desktop/Records_M2"),2);

	filteredM1 = filter_records(recordsM1, [1,2], [3,4]);
	filteredM2 = filter_records(recordsM2, [1,3], [2,4]);
	filteredC12 = filteredM1.*filteredM2;

	betasM1, betasM1_err = comp_to_paulis(filteredM1);
	betasM2, betasM2_err = comp_to_paulis(filteredM2);
	betasC12, betasC12_err = comp_to_paulis(filteredC12);
end

test_comp_to_paulis()

include("test-digitize.jl")
