using HDF5, JSON, Compat

"""
    load_data(filename)

Load hdf5 header and data
"""
function load_data(filename::AbstractString)

	h5open(filename, "r") do f

		version = read(attrs(f)["version"])[1]
		@assert version == 2

		#Some headers are malformed...
		header = try
			JSON.parse(read(f["header"])[1])
		catch
			warn("Unable to parse header JSON structure")
			Dict()
		end

		nbrDataSets = read(attrs(f)["nbrDataSets"])[1]
		data = Vector{Dict{String,Any}}(nbrDataSets)

		for ct in 1:nbrDataSets
			raw_data = read(f["DataSet$(ct)/real"]) + im*read(f["DataSet$(ct)/imag"])
			data[ct] = Dict{String,Any}("data" => raw_data)
			if has(f, "DataSet$(ct)/realvar")
				data[ct]["realvar"] = read(f["DataSet$(ct)/realvar"])
				data[ct]["imagvar"] = read(f["DataSet$(ct)/imagvar"])
				data[ct]["prodvar"] = read(f["DataSet$(ct)/prodvar"])
			end
			data[ct]["xpoints"] = read(f["DataSet$(ct)/xpoints"])
			data[ct]["xlabel"] = read(attrs(f["DataSet$(ct)/xpoints"])["label"])
			if ndims(raw_data) > 1
				data[ct]["ypoints"] = read(f["DataSet$(ct)/ypoints"])
				data[ct]["ylabel"] = read(attrs(f["DataSet$(ct)/ypoints"])["label"])
			end
			if ndims(raw_data) > 2
				data[ct]["zpoints"] = read(f["DataSet$(ct)/zpoints"])
				data[ct]["zlabel"] = read(attrs(f["DataSet$(ct)/zpoints"])["label"])
			end
		end

		# unwrap data cell array if only one data set
		if nbrDataSets == 1
			data = data[1]
		end

		return header, data
	end

end

"""
    load_data(datapath, filenum; subdir)

Search file number filenum in folder datapath/subdir. Subdir default is current date in "yymmdd"
"""
function load_data(datapath::AbstractString, filenum::Int, subdir=Dates.format(Dates.today(),"yymmdd"))
		# optionally, search for filenum instead of filename
		# search in a subdirectory with today's date, if not specified
		searchdir(path, fileid) = filter(x -> contains(x,string(filenum)) && contains(x,".h5"), readdir(path))
		filename = joinpath(datapath, subdir, searchdir(joinpath(datapath, subdir), filenum)[1])
		load_data(filename)
end
