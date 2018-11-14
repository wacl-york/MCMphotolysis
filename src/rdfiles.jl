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
  o3col = sort!(parse.(Int64, o3col))
  ifile = ["$scen.$i.txt" for i in o3col]

  # return file and folder names
  return ifile, iofolder, o3col
end #function get_fnames


"""
    getTUVdata(inpfile)

From vector `inpfile` with file names of all ozone scenarios get a dictionary
wit j value related data.
"""
function getTUVdata(inpfile)
  # Read j values
  jvals = []
  for f in inpfile
    j = readTUV(f)
    magnitude = [floor(log10(j[:jvals][i][1])) for i = 1:length(j[:jvals])]
    [j[:jvals][i] ./= 10^magnitude[i]  for i = 1:length(magnitude)]
    j[:order] = magnitude
    push!(jvals, j)
  end

  return jvals #sza, χ, TUVdata, params350, magnitude, rxns
end
