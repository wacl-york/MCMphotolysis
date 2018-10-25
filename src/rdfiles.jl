"""
    get_fnames(scen)

From the scenario name (scen) used in TUV derive names for the TUV output file
and the output folder.
"""
function setup_files(scen)
  # Test script argument and ask for input file, if missing
  if isempty(scen)
    println("Enter name for TUV scenario")
    scen = input("(Name of output file without \'.txt\'): ")
  end

  # Create output folder
  iofolder = "./params_"*strip(scen)
  try mkdir(iofolder)
  catch
    print("\033[95mFolder '$iofolder' already exists! Overwrite ")
    confirm = input("(\033[4mY\033[0m\033[95mes/\033[4mN\033[0m\033[95mo)?\033[0m ")
    if lowercase(confirm[1]) == 'y'
      cd(iofolder); files = readdir(); [rm(file) for file in files ]; cd("..")
    else println("Programme aborted. Exit code '98'."); exit(98)
    end
  end

  # Define TUV file
  ifile = scen*".txt"
  if !isfile(ifile)
    println("$ifile does not exist. Script terminated. Exit code: 99."); exit(99)
  end

  # return file and folder names
  return ifile, iofolder
end #function setup_files
