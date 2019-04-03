module MCMphotolysis

# Source directory
dir = Base.source_dir()

# Load Julia packages
using Statistics
using Dates
using ProgressMeter
using LsqFit
using PyCall
using DataFrames
using Printf
using LaTeXStrings
import LinearAlgebra.⋅
# pdf for multipage pdf plots
const pdf = PyNULL()

# Load self-made packages
import filehandling; const fh = filehandling
import filehandling.TUVdata
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


#=
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
=#


# export public functions
export j_oldpars,
       j_parameters,
       PhotData#,
       # StatData

# Import pdf from matplotlib for multipage plotting
function __init__()
  copy!(pdf, pyimport("matplotlib.backends.backend_pdf"))
end


# Include outsourced functions
include(joinpath(dir,"rdfiles.jl"))
include(joinpath(dir,"fitTUV.jl"))
include(joinpath(dir,"output.jl"))


"""
    j_oldpars(scen::String; output::Union{Bool,String}=true, DU::Number=350, MCMversion::Int64=3)

Read reactions from a `TUVdata` for a TUV output file named `scen`.txt, and store the `TUVdata`
in the following fields:

- `jval`: DataFrame with j values
- `order`: Vector with order of magnitudes
- `deg`/`:rad`: Vector with solar zenith angles in deg/rad
- `rxn`: Vector of strings with reaction labels
- `mcm`/`tuv` MCM/TUV photolysis reaction numbers
- `DU`: Ozone column from kwarg (only as info; default: `350`)


Derive MCM l,m,n-parameterisations for photolysis and return a data type `PhotData`
with the following fields:

- `l`/`m`/`n`: Vector with MCM photolysis parameters
- `sigma`: Vector of Vectors with σ-values for the 95% confidence interval for the `l`, `m`, and `n` parameters
- `converged`: Vector of booleans with `true` for a converged fit, otherwise `false`

If output is set to `true` or `"plot"` (_default_), _j_ values from TUV and the parameterisations
for all reactions are plotted to a pdf and _j_ values and errors/statistical data
are printed to a formatted text file or semi-colon separated csv file in a folder named
`params_<scen>`. The following options for output exist:

- `true` or `"plot"`: Plot TUV data with parameterisations to pdf and data to a formatted
  text file and semicolon-separated csv file
- `"data"`: Only print data to text/csv files
- `false` or `"None"`: Don't print output only return the data from the function

The `MCMversion` is needed to assign the correct MCM reaction numbers:
- `2`: MCMv3.2 or older
- `3`: MCMv3.3.1
- `4`: MCM/GECKO-A
"""
function j_oldpars(scen::String; output::Union{Bool,String}=true, DU::Number=350, MCMversion::Int64=3)
  # Initialise system time and output path/file name
  systime = now()

  # Read dataframe with j values from TUV output file
  println("load data...")
  ifile, iofolder = setup_files(scen, output)
  if ifile == ""  @info("Script stopped."); return nothing, nothing  end
  jvals = fh.readTUV(ifile, DU = DU, MCMversion = MCMversion)

  # Derive parameterisations for j values
  params = fit_jold(jvals) #, stats

  # Write output
  write_oldparams(jvals,params,iofolder,systime,output) #,stats
  plot_jold(jvals,params,systime,iofolder,output)

  return jvals, params#, stats
end #function j_oldpars

"""
    j_parameters(scen::String; output::Union{Bool,Int64,Float64,Vector{Int64},Vector{Float64}}=350, MCMversion::Int64=4)

The functions searches for a TUV output files in the current directory from the
scenario name of the TUV run `scen`. TUV files have to be of the format `<scen>.DU.txt`.

The `MCMversion` is needed to assign the correct MCM photolysis reaction numbers.
Output is printed under the following conditions, when `output` is set to:
- `false`: No output, function returns data as `TUVdata` and `PhotData`
- `true`: Output written to formatted text file `parameters.dat` and
  semicolon-separated csv file `paramters.csv`
- `Int64`/`Vector{Int64}`: Additionally to text files, parameterisations and TUV
  data are compared in pdf plots for any ozone column specified in `output`
"""
function j_parameters(scen::String;
         output::Union{Bool,Int64,Vector{Int64}}=350, MCMversion::Int64=4)
  # Initialise system time and output path/file name
  systime = now()

  # Read dataframe with j values from TUV output file
  println("load data...")
  inpfile, iofolder, O3col = getO3files(scen, output)
  if inpfile == ""  @info("Script stopped."); return nothing, nothing  end

  # Read TUV data and get original l parameters and m, n parameters for 350DU
  jvals = getTUVdata(inpfile, O3col, MCMversion)

  # Collect and fit data
  ldata, params350 = getMCMparams(jvals, O3col)
  params = fitl(ldata, jvals[1].order, O3col, params350, jvals[1].rxn)

  # Write output
  ptitle = set_titles(jvals[1])
  write_params(jvals, params, iofolder, systime, output)
  plotl(ldata, params, jvals[1].order, O3col, ptitle, iofolder, systime, output)
  plotj(jvals, params, ptitle, O3col, output, iofolder, systime)

  return jvals, params
end #function j_parameters

end # module MCMphotolysis
