using JSON, Compat, CSV, Mmap, HDF5, Dates, DataStructures, YAML

"""
    load_data(filename)

Load data file
"""
# add functionality to read new Auspex data files (format on develop branch)

function load_data(filename::AbstractString)
    datasets    = Dict{String, Array{Any}}()
    descriptors = Dict{String,Any}()
    # per qubit
    for group in readdir(filename)
        datafile = filter(x-> occursin(".dat", x), readdir(joinpath(filename, group)))
        metafile = filter(x-> occursin("meta", x), readdir(joinpath(filename, group)))
        metadata = JSON.parsefile(joinpath(filename, group, metafile[1]), dicttype=DataStructures.OrderedDict)
        open(joinpath(filename, group, datafile[1]), "r") do op
            dims = length(metadata["shape"])>1 ? tuple(metadata["shape"]...) : tuple(metadata["shape"][1],1)
            data = Mmap.mmap(op, Array{Complex{Float64},1}, prod(dims))
            data = reshape(data, reverse(dims))
            datasets[group] = data
        end
        descriptors[group] = metadata
    end
    return datasets, descriptors
end

"""
    load_data(datapath, filenum; subdir, auspex)

Search file number filenum in folder datapath/subdir. Subdir default is current date in "yymmdd"
auspex: boolean for file format (default = true)
"""
function load_data(datapath::AbstractString, filenum::Int, subdir=Dates.format(Dates.today(),"yymmdd"), h5=false)
  # optionally, search for filenum instead of filename
  # search in a subdirectory with today's date, if not specified
  #regexpr = auspex? r"(\d{4})(.h5)" : r"(\d+)(_\w+)(.h5)"
  ext = h5 ? ".h5" : ".auspex"
  regexpr = Regex("(.+)($(lpad(filenum, 4, "0")))(\\$ext)")
  files = filter(x -> occursin(regexpr, x), readdir(joinpath(datapath, subdir)))
  if length(files) == 0
    println("No matching files found!")
    return
  elseif length(files) > 1
    println("More than one matching file found!")
    return
  end
  filename = joinpath(datapath, subdir, files[1])
  h5 ? load_h5_data(filename) : load_data(filename)
end

## TO BE UPDATED
# """
# 	load_latest_data()
#     load_latest_data(ind_from_last)
#
# Load latest data file, or ind_from_last files before the latest.
# """
# function load_latest_data(logpath::AbstractString, ind_from_last::Int = 0)
# 	df = CSV.read(logpath, delim = "\t"; use_mmap = false)
# 	@assert ind_from_last < length(df.columns[1]) "Not enough data files in the history."
# 	filename = df.columns[1][end-ind_from_last]
# 	load_auspex_data(filename)
# end

"""
    load_h5_data(filename)

Load hdf5 data and descriptor in legacy auspex format

FILE STRUCTURE:
file.hdf5:
   header/ (with yml settings, if saved)
   group1/
      data/
         axis_1 (contains all tuples visited on this axis)
         axis_2 (contains all tuples visited on this axis)
         quantity_1 (the measured values)
         quantity_2 (the measured values)
      descriptor (contains references to the axes datasets below)
      axis_1 (contains the values of the axis  see note #1)
      axis_2 (contains the values of the axis - see note #1)
   group2/
      ... same structure as above ...

Note #1: Any axes referenced by the descriptor can themselves contain
         references to multiple subaxes. This indicates that the axis
         is "unstructured", i.e. that multiple coordinates where changed
         simultaneously.
"""
function load_h5_data(filename::AbstractString)


    datasets    = Dict{String, Dict{String, Vector{Any}}}()
    descriptors = Dict{String, Vector{Dict{String,Any}}}()
	header =   Dict{String, Any}()
	#info on file, such as name, ...?

	header["filename"] = filename
    h5open(filename, "r") do f
        # Find all of the group names, which will correspond to qubits when
        # using exp_factory inside of auspex. The default group name is simply
        # "main".
        group_names = names(f)
		#load measurement/instrument settings
		if "header" in group_names
			header["settings"] = YAML.load(read(attrs(f["header"])["settings"]))
		end

        for group_name in filter!(name -> name != "header", group_names)
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
    datasets, descriptors, header
end

convert_dataset(ds) = ds
# heuristic to read complex numbers stored as HDF5Compounds
function convert_dataset(ds::Vector{NamedTuple{(:r, :i), Tuple{Float32, Float32}}})
  return [complex(d.data[1], d.data[2]) for d in ds]
end
