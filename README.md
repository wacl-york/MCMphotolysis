MCMphotolysis
=============

A Julia package to retrieve updated MCM photolysis parameterisations including
a dependence on the overlying ozone column from TUV output files.


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
before your first run.


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
saved in the format `<scen>.<O3col>.txt`, where `<scen>` is a scenario name
that can be chosen freely and `<O3col` is the ozone column in DU. As TUV
allows only 6 characters for file names, `<scen>` must not be longer than
2 characters.
  
### The original MCM photolysis parameterisation

Function `j_oldpars` derives parameters for every reaction with non-zero
_j_ values in the TUV output file for the original MCM parameterisation:

    j(χ) = l·cosᵐ(χ)·exp(-n·sec(χ))


Run:

```julia
jvals, params, stats = j_oldpars("<scen>", output = true)
```

`<scen>` is the file name without `.txt`. The function creates a folder
`params_<scen>` it which it stores a text file with parameters and statistical
data for each reaction and a pdf where fits are compared to the TUV data.

By default `output` is set to true and does not have to be specified in
the above command. If set to `false`, the output folder is not created
and the function only returns the data below. Sigma is the standard error
determined with the `LsqFit` package (https://github.com/JuliaNLSolvers/LsqFit.jl.git).
Additionally, the root-mean-square error RMSE and the correlation coefficient
R² are determined and returned as `StatData`.

`jvals` of type `TUVdata` with fields:
- `jval::DataFrame`
- `order::Vector{Int64}`
- `deg::Vector{Float64}`
- `rad::Vector{Float64}`
- `rxn::Vector{String}`
- `mcm::Vector{Int64}`
- `tuv::Vector{Int64}`
- `O3col::Number`

`params` of type `PhotData` with fields:
- `l::Union{Vector{Float64},Vector{Any}}`
- `m::Vector{Float64}`
- `n::Vector{Float64}`
- `sigma::Vector{Vector{Float64}}`
- `converged::Vector{Bool}`

`stats` of type `StatData` with fields:
- `RMSE::Vector{Float64}`
- `R2::Vector{Float64}`


### Updated MCM photolysis parameterisation with ozone column dependence

Function `j_parameters` returns a refined parameterisation, which includes a dependence on the overlying ozone column. In the original parameterisation, 
`l` has been redefined to include a second order exponential dependence on the
ozone column. At 350DU, both parameterisations are identical.

    j(χ) = (l_a0 + l_b0·exp(-l_b1/[O3col])+ l_c0·exp(-l_c1/[O3col]))·cosᵐ(χ)·exp(-n·sec(χ))


Run:

```julia
jvals, params = j_parameters("M4")
```

The function works as above and by default creates a folder `params_<scen>`,
where the ozone column dependence of `l` is plotted in `lpar.pdf` and 
_j_ value fits and TUV data are compared in `jvalues.pdf`. A formatted files
`parameters.dat` lists all parameters and standard errors. Additionally,
`parameters.csv` lists the data in a csv file.

`j_parameters` does not calculate the additional statistical parameters RMSE
and R<sup>2</sup>. The parameterisation is derived by deriving the original MCM
parameterisation at 350DU, fixing `m` and `n` and recalculating `l` for every
ozone column. If values at 350DU do not exist, the next smallest ozone column
is chosen. In a next step the second order exponential decay is derived from
the ozone dependent l values.

Errors for the `l` parameters are composed of the standard errors for the
original `l` fit and the ozone column fit:

    Δlₓ = (σ(l)/l + σ(lₓ)/lₓ)·lₓ

The keyword argument `output` is slightly redefined. If set to `false`, no
output folder is created and only `TUVdata` and `PhotData` is returned.
If set to `true`, only text files are generated. To generate plots, specify
an integer or vector of integers, for which _j_ value plots are desired.
In that case, the ozone column dependence of the `l` parameter is plotted
as well.


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
