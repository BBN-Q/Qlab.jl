function read_records(baseName::String)
	#Helper function to read back in saved records

	return read_file(string(baseName, ".real")) + 1im*read_file(string(baseName, ".imag"))
end

function read_file(fileName::String)

	#Get the size of the file and open
	fileSize = filesize(fileName)
	FID = open(fileName, "r")

	#Read out the first three dimensions (recordLength x numWaveforms x numSegments)
	sizes = read(FID, Int32, 3)

	#Infer the number or records from the file size
	numRecords = int64((fileSize-13*sizeof(Int32))/prod(sizes)/sizeof(Float32))

	return mmap_array(Float32, tuple([int64(sizes), numRecords]...), FID)

end