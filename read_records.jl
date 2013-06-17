function read_records(fileName::String)
	#Helper function to read back in saved records

	return open(read_file, string(fileName, ".real")) + 1im*open(read_file, string(fileName, ".imag"))
end

function read_file(FID::IOStream)
	sizes = read(FID, Int32, 3)
	out = memio()
	ct = 0
	while !eof(FID)
		a = read(FID, Float32)
        write(out, a)
        ct += 1
    end
    seekstart(out)
    data = read(out, Float32, ct)

	bufferSize = prod(sizes)
	lastDim = fld(length(data), bufferSize)
	data = reshape(data, [int64(sizes), lastDim]...)
	return data
end