"""
    get_fnames(scen, output)

From the scenario name (scen) used in TUV, derive names for the TUV output file
and the output folder.

If `output` is `true`, create a output folder `./params_<scen>`. Ask to overwrite,
if folder already exists.
"""
function setup_files(scen, output)
  # Test script argument and ask for input file, if missing
  if isempty(scen)
    println("Enter name for TUV scenario")
    scen = input("(Name of output file without \'.txt\'): ")
  end

  # Create output folder
  iofolder = "./params_"*strip(scen)
  if output  try mkdir(iofolder)
  catch
    print("\033[95mFolder '$iofolder' already exists! Overwrite ")
    confirm = input("(\033[4mY\033[0m\033[95mes/\033[4mN\033[0m\033[95mo)?\033[0m ")
    if lowercase(confirm[1]) == 'y'
      cd(iofolder); files = readdir(); [rm(file) for file in files ]; cd("..")
    else println("Programme aborted. Exit code '98'."); exit(98)
    end
  end  end

  # Define TUV file
  ifile = scen*".txt"
  ifile = test_file(ifile)

  # return file and folder names
  return ifile, iofolder
end #function setup_files


"""
    get_O3dep_files(scen, output)

From the scenario name (`scen`) used in TUV, derive names for the TUV output file
and the output folder and return them.

If output is set, create a folder `params_<scen>` and ask to overwrite, if already
existant.
"""
function get_O3dep_files(scen, output)

  # Test script argument and ask for input file, if missing
  while !any(occursin.(scen,readdir()))
    println("Enter a valid name for a TUV scenario")
    scen = input("(Name of output file without \'.txt\'): ")
  end

  # Create output folder
  iofolder = "./params_"*strip(scen)
  if output ≠ false  try mkdir(iofolder)
  catch
    print("\033[95mFolder '$iofolder' already exists! Overwrite ")
    confirm = input("(\033[4mY\033[0m\033[95mes/\033[4mN\033[0m\033[95mo)?\033[0m ")
    if lowercase(confirm[1:1]) == "y"
      cd(iofolder); files = readdir(); for file in files  rm(file)  end; cd("..")
    else println("Programme aborted. Exit code '98'."); exit(98)
    end
  end  end

  # Define TUV file
  filelist=readdir()
  idx=findall(startswith.(filelist,scen) .& endswith.(filelist,".txt"))
  o3col=filelist[idx]
  o3col=[replace(o, "$scen." => "") for o in o3col]
  o3col=[replace(o, ".txt" => "") for o in o3col]
  o3col = sort!(parse.(Float64, o3col))
  ifile = ["$scen.$i.txt" for i in convert.(Int64, o3col)]

  # return file and folder names
  return ifile, iofolder, o3col
end #function get_fnames


function collect_TUVdata(inpfile)
  # Read j values
  jvals = []; magnitude = []; sza = Float64[]; χ = Float64[]
  params = []
  # for f in inpfile
  #   j, sza, χ = read_j(f)
  #   push!(TUVdata, j)
  # end
  for (n,f) in enumerate(inpfile)
    j = readTUV(f)
    magnitude = [floor(log10(j[:jvals][i][1])) for i = 1:length(j[:jvals])]
    for i = 1:length(magnitude)  j[:jvals][i] ./= 10^magnitude[i]  end
    j[:order] = magnitude
    push!(jvals,j)
  end
  # TUVdata = zeros(Array{Float64,3}(length(jvals[1][1]),length(jvals[1]),length(jvals)))
  # magnitude = [floor(log10(jvals[1][1,d])) for d in 1:length(jvals[1][1,:])]
  # rxns = names(jvals[1])
  # for DU = 1:length(jvals), rxn = 1:length(jvals[1])
  #   TUVdata[:,rxn,DU] = jvals[DU][rxn]/10^magnitude[rxn]
  # end
  # params350 = Matrix{Float64}(0,4)
  # j350 = DataFrame()
  # [j350[rxns[i]] = TUVdata[:,i,iDU[3]] for i = 1:length(rxns)]
  # fit, sigma, stats = fit_jold(j350,χ)
  # [params350 = vcat(params350, [fit[i].param[2] fit[i].param[3] sigma[i][2] sigma[i][3]])
  #   for i = 1:length(rxns)]

  return jvals #sza, χ, TUVdata, params350, magnitude, rxns
end
