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
Pkg.instantiate()
Pkg.precomile()
```

or go to the package manager typing `]` in the julia prompt and in the 
package manager the following commands:

```
add https://github.com/pb866/MCMphotolysis.git
instantiate
precompile
```


Usage
-----

After installation run the following lines in the REPL or write a julia
script with the following lines and call it from the terminal.

```julia
using MCMphotolysis
jvals = j_oldpars("<scenario name>")
```

If you haven't installed the package directly to your Julia default environment,
you need to activate the environment first.

Call function `j_oldpars` handing over the scenario name, i.e. the name of
the TUV input file without the `.txt` file ending. If you only want the data
returned in the REPL without file output, set the optional keyword argument
`output = false`.  
Function `j_oldpars` returns a dictionary with entries for a DataFrame with
the _j_ values (`:jvals`), where the `:order` of magnitude is stored in a
different array, arrays with the solor zenith angles (`:deg`/`:rad`),
an array with LsqFit output, statistical data (`:σ`, `:RMSE`, and `:R2`).
Moreover, a folder `params_<scenario name>` is created, where parameters are
printed to `parameters.dat` and parameterisations are visualised in 
`<scenario name>.pdf`.
 

Trouble shooting of common errors
---------------------------------

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


Version history
===============

Version 0.1.0
-------------

- Function `j_oldpars` to write the legacy MCM parameterisations and statistical
  data to a text file and plot the parameterisations in a pdf compiled in the
  folder `params_<scenario name>`
