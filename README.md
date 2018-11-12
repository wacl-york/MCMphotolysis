MCMphotolysis
=============

A Julia package to retrieve updated MCM photolysis parameterisations including
a dependence on the overlying ozone column from TUV output files.


Installation
------------

Install the Package by adding it to your environment using `Pkg`

```julia
using Pkg
Pkg.add("https://github.com/pb866/MCMphotolysis.git")
Pkg.precomile()
```

or go to the package manager typing `]` in the julia prompt and in the 
package manager the following commands:

```julia
add https://github.com/pb866/MCMphotolysis.git
precompile
```

### Trouble shooting of common errors

If you get an error message, follow the the instructions of the error message, e.g.

```julia
Pkg.build("CodecZlib")
```

If PyPlot crashes, try running Julia with the system python rather than 
the miniconda python version by rebuilding python with:

```julia
using Pkg
ENV["PYTHON"] = "path/to/python"
Pkg.build("PyPlot")
```

You can get the system python version by typing `which python` in the terminal
and copying the output.
