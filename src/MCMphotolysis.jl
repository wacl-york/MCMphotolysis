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
# PyCall python imports
# @pyimport matplotlib.backends.backend_pdf as pdf
const pdf = PyNULL()

# Load self-made packages
using filehandling
import pyp

# export public functions
export j_oldpars

# Include outsourced functions
include(joinpath(dir,"rdfiles.jl"))
include(joinpath(dir,"fitTUV.jl"))
include(joinpath(dir,"output.jl"))

function __init__()
  copy!(pdf, pyimport("matplotlib.backends.backend_pdf"))
end

"""
    j_oldpars(scen)

Read reactions from a TUV output file named `scen`.txt, derive MCM l,m,n-parameterisations
for photolysis, plot the _j_ values from TUV and the parameterisation for all reactions
to a pdf and print the j values and errors/statistical data to a text file in a folder
`params_<scen>`.
Return solar zenith angle as tuple of arrays with deg/rad, the TUV _j_ value data,
and the fit data (as array with entries for each reaction).
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
  fit, sigma, rmse, R2 = fit_jold(jvals[:jvals],jvals[:rad])
  jvals[:fit] = fit; jvals[:σ] = sigma, jvals[:RMSE] = rmse; jvals[:R2] = R2

  # # Generate output
  if output
    # println("ϑ")
    plot_jold(jvals,systime,iofolder,scen)
    # wrt_params(names(jvals),fit,magnitude,sigma,rmse,R2,iofolder,time)
  end
  return jvals
end #function j_parameters

end # module MCMphotolysis
