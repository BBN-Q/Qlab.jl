using HDF5, JSON, Compat

"""
    load_data(filename)

Load hdf5 header and data
"""
function load_data(filename::AbstractString)

	h5open(filename, "r") do f

		version = read(attrs(f)["version"])[1]
		@assert version == 2

		# Some headers are malformed...
		header = try
			JSON.parse(read(f["header"])[1])
		catch
			warn("Unable to parse header JSON structure")
			Dict{String,Any}()
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
    load_auspex_data(filename)

Load hdf5 data and descriptor in auspex format

FILE STRUCTURE:
file.hdf5:
   group1/
      data/
         axis_1 (contains all tuples visited on this axis)
         axis_2 (contains all tuples visited on this axis)
         quantity_1 (the measured values)
         quantity_2 (the measured values)
      descriptor (contains references to the axes datasets below)
      axis_1 (contains the values of the axis — see note #1)
      axis_2 (contains the values of the axis - see note #1)
   group2/
      ... same structure as above ...

Note #1: Any axes referenced by the descriptor can themselves contain
         references to multiple subaxes. This indicates that the axis
         is "unstructured", i.e. that multiple coordinates where changed
         simultaneously.
"""
function load_auspex_data(filename::AbstractString)


    datasets    = Dict{String, Dict{String, Vector{Any}}}()
    descriptors = Dict{String, Vector{Dict{String,Any}}}()

    h5open(filename, "r") do f
        
        # Find all of the group names, which will correspond to qubits when
        # using exp_factory inside of auspex. The default group name is simply
        # "main".
        group_names = names(f)

        for group_name in group_names
            g = f[group_name]

            # Read in the descriptor
            desc_refs = read(g["descriptor"])

            desc = Vector{Dict{String,Any}}()
            for desc_ref in desc_refs
                axis = Dict{String,Any}()
                if reinterpret(Bool, read(attrs(g[desc_ref])["unstructured"]))
                    axis_refs      = read(g[desc_ref])
                    axis["points"] = [read(g[axis_ref]) for axis_ref in axis_refs]
                    axis["name"]   = [read(attrs(g[axis_ref])["name"]) for axis_ref in axis_refs]
                    axis["unit"]   = [read(attrs(g[axis_ref])["unit"]) for axis_ref in axis_refs]
                else
                    axis["points"] = read(g[desc_ref])
                    axis["name"]   = read(attrs(g[desc_ref])["name"])
                    axis["unit"]   = read(attrs(g[desc_ref])["unit"])
                end
                # See if we need to get metadata:
                if "metadata" in names(attrs(g[desc_ref]))
                    meta_ref = read(attrs(g[desc_ref])["metadata"])
                    meta_enum_ref = read(attrs(g[desc_ref])["metadata_enum"])
                    axis["metadata"] = [read(g[meta_enum_ref])[m+1] for m in read(g[meta_ref])]
                end
                
                push!(desc, axis)
            end

            descriptors[group_name] = desc

            # Read in the actual data
            datacol_names = names(g["data"])
            datasets[group_name] = Dict{String, Vector{Any}}()
            for datacol_name in datacol_names
              dataset = read(g["data"][datacol_name])
              # heuristic to read complex numbers stored as HDF5Compounds
              datasets[group_name][datacol_name] = convert_dataset(dataset)
            end
        end
    end
    datasets, descriptors
end

convert_dataset(ds) = ds
# heuristic to read complex numbers stored as HDF5Compounds
function convert_dataset(ds::Vector{HDF5.HDF5Compound{2}})
  return [complex(d.data[1], d.data[2]) for d in ds]
end

"""
    load_data(datapath, filenum; subdir)

Search file number filenum in folder datapath/subdir. Subdir default is current date in "yymmdd"
"""
function load_data(datapath::AbstractString, filenum::Int, subdir=Dates.format(Dates.today(),"yymmdd"), auspex::Bool=true)
  # optionally, search for filenum instead of filename
  # search in a subdirectory with today's date, if not specified
  filename = joinpath(datapath, subdir, filter(x -> contains(x,string(filenum)) && contains(x,".h5"), readdir(joinpath(datapath, subdir)))[1])
  auspex? load_auspex_data(filename) : load_data(filename)
end

"""
    load_auspex_data(datapath, filenum; subdir)

Search file number filenum in folder datapath/subdir. Subdir default is current date in "yymmdd"
"""
function load_auspex_data(datapath::AbstractString, filenum::Int, subdir=Dates.format(Dates.today(),"yymmdd"))
	# optionally, search for filenum instead of filename
	# search in a subdirectory with today's date, if not specified
	filename = joinpath(datapath, subdir, filter(x -> contains(x,string(filenum)) && contains(x,".h5"), readdir(joinpath(datapath, subdir)))[1])
	load_auspex_data(filename)
end
