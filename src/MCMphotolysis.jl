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
using PyPlot, PyCall
using DataFrames
using Printf
# PyCall python imports
# @pyimport matplotlib.backends.backend_pdf as pdf
const pdf = PyNULL()

# Load self-made packages
using filehandling
import pyp

# export public functions
export j_oldpars,
       j_parameters

# Include outsourced functions
include(joinpath(dir,"rdfiles.jl"))
include(joinpath(dir,"fitTUV.jl"))
include(joinpath(dir,"output.jl"))

function __init__()
  copy!(pdf, pyimport("matplotlib.backends.backend_pdf"))
end

"""
    j_oldpars(scen::String; output::Bool=true)

Read reactions from a TUV output file named `scen`.txt, derive MCM l,m,n-parameterisations
for photolysis and return a dictionary with the following entries:

- `:jvals`: DataFrame with j values (devided by maximum order of magnitude, i.e. without the `e-...`)
- `:order`: Vector with order of magnitudes
- `:deg`/`:rad`: Vector with solar zenith angles in deg/rad
- `:fit`: Vector with LsqFit data
- `:σ`: Vector of Vectors with σ-values for the 95% confidence interval for the `l`, `m`, and `n` parameters
- `:RMSE`: Vector of root mean square errors
- `:R2`: Vector of correlation cofficients R²

If output is set to `true` (_default_), _j_ values from TUV and the parameterisations
for all reactions are plotted to a pdf and _j_ values and errors/statistical data
are printed to a text file in a folder named `params_<scen>`.
"""
function j_oldpars(scen::String; output::Bool=true)
  # Initialise system time and output path/file name
  systime = now()

  # Read dataframe with j values from TUV output file
  println("load data...")
  ifile, iofolder = setup_files(scen, output)
  jvals = readTUV(ifile)

  magnitude = [floor(log10(jvals[:jvals][i][1])) for i = 1:length(jvals[:jvals])]
  for i = 1:length(magnitude)  jvals[:jvals][i] ./= 10^magnitude[i]  end
  jvals[:order] = magnitude

  # Derive parameterisations for j values
  fitstats = fit_jold(jvals[:jvals],jvals[:rad])
  jvals = merge(jvals, fitstats)

  # # Generate output
  if output
    # println("ϑ")
    plot_jold(jvals,systime,iofolder,scen)
    wrt_params(jvals,iofolder,systime)
  end
  return jvals
end #function j_parameters


"""
    j_parameters(scen::String; output::Union{Bool,Int64,Float64,Vector{Int64},Vector{Float64}}=350)

The functions searches for a TUV output file in the current directory from the scenario name of
the TUV run `scen` (output file name without `.txt`) and create a folder `params_<scen>` with a
file `parameters.dat` listing the fitting parameters for the MCM photolysis parameterisations and
`<scen>.pdf` with a graphical display of the TUV calculated data and the MCM parameterisation.
"""
function j_parameters(scen::String;
         output::Union{Bool,Int64,Float64,Vector{Int64},Vector{Float64}}=350)
  # Initialise system time and output path/file name
  systime = now()

  # Read dataframe with j values from TUV output file
  println("load data...")
  inpfile, iofolder, o3col = get_O3dep_files(scen, output)

  # Read TUV data and get original l parameters and m, n parameters for 350DU
  # Data is rescaled to exclude the order of magnitude
  # sza, χ, TUVdata, params350, magnitude, rxns =
  jvals = collect_TUVdata(inpfile)

  # parMCM, sigMCM, jMCM, fit = fit_j(TUVdata, params350, o3col, χ, iDU, rxns)
  #
  # plot_j(sza,χ,o3col,TUVdata,jMCM,magnitude,rxns,iDU,time,iofolder,scen)

  return jvals #sza, TUVdata, fit
end #function j_parameters

end # module MCMphotolysis
