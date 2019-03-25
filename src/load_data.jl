using JSON, Compat, CSV, Mmap

"""
    load_data(filename)

Load data file
"""
# add functionality to read new Auspex data files (format on develop branch)

function load_data(filename::AbstractString)
    datasets    = Dict{String, Matrix{Any}}()
    descriptors = Dict{String,Any}()
    # per qubit
    for group in readdir(filename)
        datafile = filter(x-> occursin(".dat", x), readdir(joinpath(filename, group)))
        metafile = filter(x-> occursin("meta", x), readdir(joinpath(filename, group)))
        metadata = JSON.parsefile(joinpath(filename, group, metafile[1]))
        open(joinpath(filename, group, datafile[1]), "r") do op
            if length(metadata["shape"])>1
                dims = tuple(metadata["shape"]...)
            else
                dims = tuple(metadata["shape"][1],1)
            end
            data = Mmap.mmap(op, Matrix{Complex{Float64}}, dims)
            datasets[group] = data
        end
        descriptors[group] = metadata
    end
    return datasets, descriptors
end

###TO BE UPDATED
"""
    load_data(datapath, filenum; subdir, auspex)

Search file number filenum in folder datapath/subdir. Subdir default is current date in "yymmdd"
auspex: boolean for file format (default = true)
"""
function load_data(datapath::AbstractString, filenum::Int, subdir=Dates.format(Dates.today(),"yymmdd"), auspex=true)
  # optionally, search for filenum instead of filename
  # search in a subdirectory with today's date, if not specified
  #regexpr = auspex? r"(\d{4})(.h5)" : r"(\d+)(_\w+)(.h5)"
  regexpr = auspex ? Regex("(.+)($(lpad(filenum, 4, "0")))(\\.h5)") : Regex("($filenum)(\\_.+)(\\.h5)")
  files = filter(x -> ismatch(regexpr, x), readdir(joinpath(datapath, subdir)))
  if length(files) == 0
    println("No matching files found!")
    return
  elseif length(files) > 1
    println("More than one matching file found!")
    return
  end
  filename = joinpath(datapath, subdir, files[1])
  auspex ? load_auspex_data(filename) : load_data(filename)
end

"""
	load_latest_data()
    load_latest_data(ind_from_last)

Load latest data file, or ind_from_last files before the latest.
"""
function load_latest_data(logpath::AbstractString, ind_from_last::Int = 0)
	df = CSV.read(logpath, delim = "\t"; use_mmap = false)
	@assert ind_from_last < length(df.columns[1]) "Not enough data files in the history."
	filename = df.columns[1][end-ind_from_last]
	load_auspex_data(filename)
end
