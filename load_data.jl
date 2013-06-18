using HDF5, JSON

function load_data(filename)

	h5open(filename, "r") do f

		version = read(attrs(f)["version"])[1]
		@assert version == 2

		header = parse(read(f["header"])[1])
		nbrDataSets = read(attrs(f)["nbrDataSets"])[1]
		data = cell(nbrDataSets)

		for ct in 1:nbrDataSets
			raw_data = read(f["DataSet$(ct)/real"]) + im*read(f["DataSet$(ct)/real"])
			data[ct] = {"data" => raw_data}
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

		header, data
	end

end