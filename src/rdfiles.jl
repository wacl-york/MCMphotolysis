"""
    get_fnames(scen, output)

From the scenario name (scen) used in TUV, derive names for the TUV output file
and the output folder.

If `output` is `true`, create a output folder `./params_<scen>`. Ask to overwrite,
if folder already exists.
"""
function setup_files(scen, output)
  # Test for existence of input file
  ifile = scen*".txt"
  while !isfile(ifile)
    println("Wrong or missing TUV scenario name.")
    print("(Enter name of output file without \'.txt\') or press <ENTER> to quit: ")
    scen = readline()
    if scen == ""  return "", nothing  end
    ifile = scen*".txt"
  end

  # Create output folder
  iofolder = "./params_"*strip(scen)
  if output ≠ false && output ≠ "None"  try mkdir(iofolder)
  catch
    print(string("\033[93m\33[1mWarning:\33[0m\033[93m Folder '$iofolder'",
      " already exists! Overwrite ",
      "(\033[4mY\033[0m\033[93mes/\033[4mN\033[0m\033[93mo)?\033[0m "))
    confirm = readline()
    if !isempty(confirm) && lowercase(confirm[1]) == 'y'
      cd(iofolder); files = readdir(); [rm(file) for file in files ]; cd("..")
    else println("Stop function 'setup_files'."); return "", nothing
    end
  end  end

  # return file and folder names
  return ifile, iofolder
end #function setup_files


"""
    getO3files(scen, output)

From the scenario name (`scen`) used in TUV, derive names for the TUV output file
and the output folder and return them.

If output is set, create a folder `params_<scen>` and ask to overwrite, if already
existant.
"""
function getO3files(scen, output)

  # Test script argument and ask for input file, if missing
  while !any(occursin.(scen,readdir()))
    println("Enter a valid name for a TUV scenario")
    print("(Name of output file without \'.txt\' or stop by pressing <ENTER>): ")
    scen = readline()
  end

  # Create output folder
  iofolder = "./params_"*strip(scen)
  if output ≠ false  try mkdir(iofolder)
  catch
    print(string("\033[93m\33[1mWarning:\33[0m\033[93m Folder '$iofolder'",
      " already exists! Overwrite ",
      "(\033[4mY\033[0m\033[93mes/\033[4mN\033[0m\033[93mo)?\033[0m "))
    confirm = readline()
    if !isempty(confirm) && lowercase(confirm[1:1]) == "y"
      cd(iofolder); files = readdir(); for file in files  rm(file)  end; cd("..")
    else println("Stop function 'getO3files'."); return "", nothing, nothing
    end
  end  end

  # Define TUV file
  filelist=readdir()
  idx=findall(startswith.(filelist,scen*".") .& endswith.(filelist,".txt"))
  o3col=filelist[idx]
  o3col=[replace(o, "$scen." => "") for o in o3col]
  o3col=[replace(o, ".txt" => "") for o in o3col]
  o3col = sort!(parse.(Int64, o3col))
  ifile = ["$scen.$i.txt" for i in o3col]
  if isempty(ifile)  ifile = ""  end

  # return file and folder names
  return ifile, iofolder, o3col
end #function get_fnames


"""
    getTUVdata(inpfile, O3col, MCMversion)

From vector `inpfile` with file names of all ozone scenarios and vector `O3col`
with all ozone column values, get a vector of `TUVdata` with j value related data.
"""
function getTUVdata(inpfile, O3col, MCMversion)
  # Read j values
  jvals = []
  for i = 1:length(inpfile)
    j = fh.readTUV(inpfile[i], DU = O3col[i], MCMversion = MCMversion)
    push!(jvals, j)
  end

  return jvals #sza, χ, TUVdata, params350, magnitude, rxns
end
