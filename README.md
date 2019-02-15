MCMphotolysis
=============

A Julia package to retrieve updated MCM photolysis parameterisations including
a dependence on the overlying ozone column from TUV output files.
The package is tested against MacOS 10.14.1.


Installation
------------

Install the Package by adding it to your environment going to the package
manager typing `]` in the julia prompt and in the package manager the following commands:

```julia
julia> ]
pkg> add https://github.com/pb866/filehandling.git
pkg> add https://github.com/pb866/pyp.git
pkg> add https://github.com/pb866/MCMphotolysis.git
pkg> build PyPlot
pkg> instantiate
pkg> precompile
```

`MCMphotolysis` relies on two unregistered packages, which have to be installed
beforehand. If you haven't installed PyPlot previously, you should build it
before your first run or set plot output to false.


Usage
-----

After installation, go back to the julia prompt by hitting backspace and
import the package with `using MCMphotolysis`.
If you haven't installed the package directly to your Julia default environment,
you need to activate the environment first.

The package has two functions:
- `j_oldpars` for fitting the original l,m,n MCM photolysis parameterisations
- `j_parameters` for fitting updated l<sub>a0</sub>,l<sub>b0</sub>,l<sub>b1</sub>,
  l<sub>c0</sub>,l<sub>c1</sub>,m,n MCM photolysis parameterisations with a
  dependence on the overlying ozone column.

Both functions rely on TUV 5.2 (or version with same format) output files,
with j for every reaction at the chosen height level. For the updated
parameterisations, you need TUV files for every ozone column in the chosen
range (derived for 50 to 600DU every 50DU in the MCM). Files need to be
saved in the format `<scen>.<DU>.txt`, where `<scen>` is a scenario name
that can be chosen freely and `<DU` is the ozone column in Dobson units. As TUV
allows only 6 characters for file names, `<scen>` must not be longer than
2 characters.

### Deriving MCM parameterisations

Function `j_oldpars` derives parameters for every reaction with non-zero
_j_ values in the TUV output file for the original MCM parameterisation:

    j(χ) = l·cosᵐ(χ)·exp(-n·sec(χ))

In function `j_parameters`, l as been re-defined as a function of ozone column:

    l([DU]) = l_a0 + l_b0·exp(-l_b1/[DU])+ l_c0·exp(-l_c1/[DU])


Run:

```julia
jvals, params = j_oldpars("<scen>", **kwargs)
jvals, params = j_parameters("<scen>", **kwargs)
```

`<scen>` is the file name without `.txt` and without the ozone column in j_parameters.
The function creates a folder `params_<scen>` in which it stores a text files
(formatted and semicolon-separated) with parameters and 95% sigma values for
every reaction and a pdf where fits are compared to the TUV data.

By default all `output` is saved, but it can be (partially) switched off by the
following options:

`j_oldpars`:
- `true` or `"plot"`: all output is saved
- `"data"`: only text files are saved, no pdf with plots
- `false` or `"None"`: no output, only function return values

`j_parameters`:
- `Int64`/`Vector{Int64}`: text files and pdf with plots for the ozone columns
  specified in output either as integer or vector of several integers
- `true`: only text files are saved, no pdf with plots
- `false`: no output, only function return values

In addition to the output, the function returns to data types `TUVdata` and `PhotData`
with the following fields:

`TUVdata`:
- `jval::DataFrame` with j values of all reactions
- `order::Vector{Int64}` with order of magnitude for all j values for plot formatting
- `deg::Vector{Float64}` with solar zenith angles in deg
- `rad::Vector{Float64}` with solar zenith angles in rad
- `rxn::Vector{String}` with vector of TUV reaction strings
- `mcm::Vector{Int64}` with MCM photolysis reaction numbers
- `tuv::Vector{Int64}` with TUV photolysis reaction numbers
- `O3col::Number` with ozone column value in current scenario as specified by kwarg

`PhotData`:
- `l::Union{Vector{Float64},Vector{Any}}` with MCM l parameters for every reaction
- `m::Vector{Float64}` with MCM m parameters for every reaction
- `n::Vector{Float64}` with MCM n parameters for every reaction
- `sigma::Vector{Vector{Float64}}` with 95% sigma values for every parameter in every reaction
- `converged::Vector{Bool}` with flags for convergence

Sigma is the standard error determined with the `LsqFit` package (https://github.com/JuliaNLSolvers/LsqFit.jl.git).

Errors for the updated `l` parameters in `j_parameters` are composed of the
standard errors for the original `l` fit and the ozone column fit:

    Δlₓ = (σ(l)/l + σ(lₓ)/lₓ)·lₓ


Trouble shooting of common errors
---------------------------------

If you get an error message, follow the the instructions of the error message, e.g.

```julia
Pkg.build("CodecZlib")
```

If PyPlot crashes or fails to install, try running Julia with the system python rather than
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

Version 0.3.0
-------------
- Updated database files and new kwarg `MCMversion` for correct MCM reaction numbers
- Remove `StatData`; `RMSE` and `R2` are not supported in v0.3.0
- New options to only print output to text files in `j_oldpars`; additional
  semicolon-separated output
- Rename kwarg `O3col` to `DU`
- On failure, don't exit Julia, only abort function and return nothing
- Rename internal functions
- Code clean-up and bug fixes

Version 0.2.1
-------------
- Improved error handling: assign `Inf` to sigmal, if conversion fails to avoid
  errors from LsqFit and the abortion of the script
- Fix #4

Version 0.2.0
-------------
- Additional function `j_parameters` with improved parameterisations including
  a dependence on the overlying ozone column. Additional csv file and pdf with
  the ozone column dependence fits (lpar.pdf)

Version 0.1.0
-------------
- Function `j_oldpars` to write the legacy MCM parameterisations and statistical
  data to a text file and plot the parameterisations in a pdf compiled in the
  folder `params_<scenario name>`
