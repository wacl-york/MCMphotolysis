module MCMphotolysis

# Source directory
dir = Base.source_dir()

# Load Julia packages
using Statistics
using LinearAlgebra
using Juno: input #get terminal input inside or outside Atom
using Dates
using ProgressMeter
using LsqFit
using PyCall
using DataFrames
using Printf
using LaTeXStrings
# PyCall python imports
const pdf = PyNULL()

# Load self-made packages
using filehandling
import pyp

### NEW TYPES
"""
    struct PhotData

`PhotData` has the following fields:
- `l::Union{Vector{Float64},Vector{Vector{Float64}}}`: Either a vector of floats for the
  original parameterisation or a vector of vectors of floats for the updated `l` parameter
  in the MCM photolysis parameterisation in every reaction of TUV
- `m::Vector{Float64}`: Vector of floats with `m` parameters of the MCM photolysis parameterisation
  for every reaction in TUV
- `n::Vector{Float64}`: Vector of floats with `n` parameters of the MCM photolysis parameterisation
  for every reaction in TUV
- `sigma::Vector{Vector{Float64}}`: Vector with vectors of floats holding all standard
  deviations for every parameter of the MCM photolysis parameterisation in every reaction
  of TUV in the order l or updated l_0, l_1, l_2, l_3, and l_4, m, and n
- `converged::Vector{Bool}`: Flags, whether the fit in each reaction has converged
"""
struct PhotData
  l::Union{Vector{Float64},Vector{Any}}
  m::Vector{Float64}
  n::Vector{Float64}
  sigma::Vector{Vector{Float64}}
  converged::Vector{Bool}
end


"""
    struct StatData

Immutable structure with additional statistical data `RMSE` and `R2` (correlation factor)
for the original MCM photolysis parameterisation. Data is stored in vectors of floats
for every TUV reaction.
"""
struct StatData
  RMSE::Vector{Float64}
  R2::Vector{Float64}
end


# export public functions
export j_oldpars,
       j_parameters,
       PhotData


function __init__()
  copy!(pdf, pyimport("matplotlib.backends.backend_pdf"))
end


# Include outsourced functions
include(joinpath(dir,"rdfiles.jl"))
include(joinpath(dir,"fitTUV.jl"))
include(joinpath(dir,"output.jl"))


"""
    j_oldpars(scen::String; output::Bool=true, O3col::Number=350)

Read reactions from a TUV output file named `scen`.txt, derive MCM l,m,n-parameterisations
for photolysis and return a data type `PhotData` with the following fields:

- `jval`: DataFrame with j values (devided by maximum order of magnitude, i.e. without the `e-...`)
- `order`: Vector with order of magnitudes
- `deg`/`:rad`: Vector with solar zenith angles in deg/rad
- `rxn`: Vector of strings with reaction labels
- `O3col`: Ozone column from kwarg
- `l`/`m`/`n`: Vector with MCM photolysis parameters
- `sigma`: Vector of Vectors with σ-values for the 95% confidence interval for the `l`, `m`, and `n` parameters
- `RMSE`: Vector of root mean square errors
- `R2`: Vector of correlation cofficients R²
- `converged`: Vector of booleans with `true` for a converged fit, otherwise `false`

If output is set to `true` (_default_), _j_ values from TUV and the parameterisations
for all reactions are plotted to a pdf and _j_ values and errors/statistical data
are printed to a text file in a folder named `params_<scen>`.
"""
function j_oldpars(scen::String; output::Bool=true, O3col::Number=350)
  # Initialise system time and output path/file name
  systime = now()

  # Read dataframe with j values from TUV output file
  println("load data...")
  ifile, iofolder = setup_files(scen, output)
  jvals = readTUV(ifile, O3col)

  # Derive parameterisations for j values
  params, stats = fit_jold(jvals)

  # Write output
  if output
    plot_jold(jvals,params,systime,iofolder)
    wrt_params(jvals,params,stats,iofolder,systime)
  end
  return jvals, params, stats
end #function j_oldpars

"""
    j_parameters(scen::String; output::Union{Bool,Int64,Float64,Vector{Int64},Vector{Float64}}=350)

The functions searches for a TUV output file in the current directory from the scenario name of
the TUV run `scen` (output file name without `.txt`) and create a folder `params_<scen>` with a
file `parameters.dat` listing the fitting parameters for the MCM photolysis parameterisations and
`<scen>.pdf` with a graphical display of the TUV calculated data and the MCM parameterisation.
"""
function j_parameters(scen::String;
         output::Union{Bool,Int64,Vector{Int64}}=350)
  # Initialise system time and output path/file name
  systime = now()

  # Read dataframe with j values from TUV output file
  println("load data...")
  inpfile, iofolder, O3col = getO3files(scen, output)

  # Read TUV data and get original l parameters and m, n parameters for 350DU
  jvals = getTUVdata(inpfile, O3col)

  # Collect and fit data
  ldata, params350 = getMCMparams(jvals, O3col)
  params = fitl(ldata, jvals[1].order, O3col, params350, jvals[1].rxn)

  # Write output
  ptitle = beautify_chem(jvals[1].rxn)
  wrt_newparams(jvals, params, iofolder, systime, output)
  plotl(ldata, params, jvals[1].order, O3col, ptitle, iofolder, systime, output)
  plotj(jvals, params, ptitle, O3col, output, iofolder, systime)

  return jvals, params
end #function j_parameters

end # module MCMphotolysis
