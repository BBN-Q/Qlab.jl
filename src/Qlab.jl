module Qlab
	using Compat, LaTeXStrings

	export load_data, KT_estimation, pauliAlphabet, filter_records, comp_to_paulis, savitsky_golay, read_records, digitize

	include("load_data.jl")
	include("KT_estimation.jl")
	include("comp_to_paulis.jl")
	include("SavitskyGolay.jl")
	include("read_records.jl")
	include("digitize.jl")
	include("common_fits.jl")
	include("tomography.jl")
	include("plotting.jl")
	include("shortcuts.jl")
end
