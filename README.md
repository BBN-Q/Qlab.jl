![Build Status](https://github.com/BBN-Q/Qlab.jl/workflows/CI/badge.svg)

[![codecov](https://codecov.io/gh/BBN-Q/Qlab.jl/branch/master/graph/badge.svg?token=MvxAbHBcKL)](https://codecov.io/gh/BBN-Q/Qlab.jl)

Qlab.jl
==========

Data manipulation and analysis tools tailored for quantum computing experiments in conjunction with [Auspex](https://github.com/BBN-Q/Auspex.git).  Currently working with Julia v1.0.

## Installation

```
(v1.3) pkg> add https://github.com/BBN-Q/Qlab.jl
```
The code base also uses some system tools and python libraries for building libraries and plotting data with PyPlot.jl.  You'll want to make sure your system has these.

In CentOS:
```bash
yum -y install epel-release
yum install gcc gcc-c++ make bzip2 hdf5 libaec libgfortran libquadmath
```
In Ubuntu/Debian:
```bash
apt-get install gcc g++ gcc-7-base make libaec0 libgfortran4 libhdf5-100 libquadmath0 libsz2
```

### Python

You'll need a working version of PyPlot.  In some cases the package manager has trouble getting this right on all systems/OSs.  If you run into issues, we recommend using Conda.jl manually:
```julia
using Pkg
Pkg.add("PyCall")
Pkg.add("Conda")
ENV["PYTHON"] = ""
Pkg.build("PyCall")
using Conda
Conda.add("matplotlib")
Conda.add("seaborn")
Pkg.add("PyPlot")
```
In most cases, Julia should take care of this for you.

### Other dependancies

Qlab.jl depends on several other Julia packages that have biniary dependencies.  These should mostly be taken care of by the package manager.  One important exception is HDF5 and its libhdf5 dependancy.  This library manages the handling of HDF5 files and is currently maintained for backwards compatibility.  The version of libhdf5 which produced any data files you want to analyze must match the library version used to create the files.  You may need to add the path the the right version of libhdf5 to the Libdl path in Julia and rebuild HDF5:
```julia
push!(Libdl.DL_LOAD_PATH, "/opt/local/lib")
Pkg.build("HDF5")
```
where `/opt/local/lib` is the path to the correct version of libhdf5.  See the documentation from HDF5.jl for more details.  Currently only version of hdf5 1.8.2 - 1.8.17 are supported.  If you're not planning to use HDF5 files, you shouldn't have to worry about the library versions matching.

## Copyright

Raytheon BBN Technologies.

## License

Apache Lincense 2.0 ([summary](https://tldrlegal.com/license/apache-license-2.0-(apache-2.0)))
