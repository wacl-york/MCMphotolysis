module MCMphotolysis

# Source directory
dir = Base.source_dir()

# Track changes during development
using Pkg
Pkg.activate(".")
using Revise

# Load Julia packages
using Statistics
using Juno: input #get terminal input inside or outside Atom
using Dates: now
using ProgressMeter
using LsqFit

# Load self-made packages
using filehandling

# export public functions
export j_oldpars

# Include outsourced functions
include(joinpath(dir,"rdfiles.jl"))
include(joinpath(dir,"fitTUV.jl"))

"""
    j_oldpars(scen)

Read reactions from a TUV output file named `scen`.txt, derive MCM l,m,n-parameterisations
for photolysis, plot the _j_ values from TUV and the parameterisation for all reactions
to a pdf and print the j values and errors/statistical data to a text file in a folder
`params_<scen>`.
Return solar zenith angle as tuple of arrays with deg/rad, the TUV _j_ value data,
and the fit data (as array with entries for each reaction).
"""
function j_oldpars(scen)
  # Initialise system time and output path/file name
  systime = now()

  # Read dataframe with j values from TUV output file
  println("load data...")
  ifile, iofolder = setup_files(scen)
  jvals, sza, χ = readTUV(ifile)

  magnitude = [floor(log10(jvals[i][1])) for i = 1:length(jvals)]
  for i = 1:length(magnitude)  jvals[i] ./= 10^magnitude[i]  end

  # Derive parameterisations for j values
  fit, sigma, rmse, R2 = fit_jold(jvals,χ)

  #=
  # Generate output
  plot_jold(sza,χ,jvals,magnitude,fit,time,iofolder,scen)
  wrt_params(names(jvals),fit,magnitude,sigma,rmse,R2,iofolder,time)

  # Write data to output text file

  =#
  return (sza, χ), jvals, fit, (sigma, rmse, R2)
end #function j_parameters

end # module MCMphotolysis

#=
systime = Dates.now()

println("load data...")
scen = "testM4"
ifile, iofolder = setup_files(scen)
jvals, sza, χ = readTUV(ifile)

scen = Juno.input("what? ")
=#
