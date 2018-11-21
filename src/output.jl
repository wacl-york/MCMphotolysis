"""
    plot_j(jvals::Dict{Symbol,Any},systime::DateTime,iofolder::String,ifile::String)

Plot j values saved in `jvals[:jvals]` and parameterisations derived from parameters
in `jvals[:fit]` to file `ifile.pdf` in the `iofolder` together with the
`systime` of creation.
"""
function plot_jold(jvals::filehandling.TUVdata, params::PhotData,systime::DateTime,
  iofolder::String,ifile::String)

  # Format titles with Latex using header of jvals
  ptitle = beautify_chem(jvals.rxn)

  # Initialise array of plots and define x data
  opfile = pdf[:PdfPages]("$iofolder/$ifile.pdf")

  # Loop over all reactions
  @showprogress 1 "plot data..." for i=1:length(params.l)
    # define parameters
    # Calculate parameterised values
    χ = collect(0:π/360:π/2)
    jpar = params.l[i]./10^jvals.order[i].⋅cos.(χ).^params.m[i].⋅exp.(-params.n[i].⋅sec.(χ))
    # Load TUV data and parameterisation for plotting
    jplt = pyp.load_PlotData(DataFrame(χ=jvals.deg, j=jvals.jval[i]./10^jvals.order[i]),
                            label = "TUV data", pt = "s", lc = "black", lt = "None")
    pplt = pyp.load_PlotData(DataFrame(χ=collect(0:0.5:90), j=jpar), lc = "red", lw = 2,
                            label="MCM Parameterisaton")
    # Plot TUV data

    fig, ax = pyp.plot_data(jplt, pplt, ti = ptitle[i],
    xlabel = "solar zenith angle χ",
    ylabel = "j / 10\$^{$(Int(jvals.order[i]))}\$ s\$^{-1}\$",
    xlims=(0,90), ylims=(0,nothing), maj_xticks=15, min_xticks=5, ti_offset=-2,
    legpos="lower left")
    # Plot time stamp
    ax[:annotate]("created $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))",
          xy=[0.02;0.22], xycoords="axes fraction", ha="left", va="bottom")

    # save plot to temporary png
    opfile[:savefig](fig)
    close()
  end

  # compile pngs in single pdf and delete pngs
  opfile[:close]()
end #function plot_j


"""


"""
function plotl(ldat, params, order, O3col, ptitle, iofolder, systime, output)

  # Only print, if output is set to true
  if output == false  return  end

  # Initialise array of plots and define x data
  opfile = pdf[:PdfPages]("$iofolder/lpar.pdf")
  # Loop over reactions
  @showprogress 1 "plot l..." for i = 1:length(ptitle)
    # define parameters
    # Calculate parameterised values
    o3 = collect(O3col[1]:O3col[end])
    lpar = lnew(o3, params.l[i])/10^order[i] #lnew is defined in fitTUV.jl
    # Load TUV data and parameterisation for plotting
    ldata  = pyp.load_PlotData(DataFrame(o3=O3col, l=ldat[i,:]./10^order[i]),
                 label = "l data", pt = "s", lc = "black", lt = "None")
    lparam = pyp.load_PlotData(DataFrame(o3=o3, l=lpar), lc = "red", lw = 2,
                 label="fit")
    # Plot TUV data

    fig, ax = pyp.plot_data(ldata, lparam, ti = ptitle[i],
    xlabel = "ozone column / DU",
    ylabel = "l / 10\$^{$(Int(order[i]))}\\,\$s\$^{-1}\$",
    xlims=(O3col[1],O3col[end]),ti_offset=-2,
    legpos="upper right")
    # Plot time stamp
    ax[:annotate]("created $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))",
          xy=[0.2;0.9], xycoords="axes fraction", ha="left", va="bottom")

    # save plot to temporary png
    opfile[:savefig](fig)
    close()
  end

  # compile pngs in single pdf and delete pngs
  opfile[:close]()
end


"""
    wrt_params(jvals, iofolder, systime)

Write the parameters, statistical data, and reaction labels stored in the
dictionary `jvals` to the file `parameters.dat` in the `iofolder` and state
the `systime` of creation.
"""
function wrt_params(jvals, params, stats, iofolder, systime)

  # Open output file
  open("$iofolder/parameters.dat","w") do f
    # Print header
    println(f,
    "Parameters and statistical data for parameterisation of photolysis processes")
    println(f,
    "in the Master Chemical Mechanism (MCM; http://mcm.york.ac.uk/) using:")
    println(f, "\nj / s-1 = l·(cos(x))^m·exp(-n·sec(x))\n")
    println(f, "Fits and standard errors are derived with Julia package LsqFit v0.6.0")
    println(f, "(https://github.com/JuliaNLSolvers/LsqFit.jl.git).")
    println(f, "created $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))")
    println(f,"\n                 P a r a m e t e r s               S t a t i s t i c s")
    println(f,"     l / s-1              m              n         RMSE / s-1    R^2      Reaction")

    # Loop over reactions
    for i = 1:length(jvals.rxn)
      # Print parameters, statistical data, and reaction label to output file
      @printf(f,"(%6.3f±%.3f)e%d    %.3f±%.3f    %.3f±%.3f    %.3e    %.4f    %s\n",
      params.l[i]/10^jvals.order[i], params.sigma[i][1]/10^jvals.order[i], Int(jvals.order[i]),
      params.m[i], params.sigma[i][2], params.n[i], params.sigma[i][3],
      stats.RMSE[i], stats.R2[i], jvals.rxn[i])
    end
  end
end # function wrt_params


function wrt_newparams(jvals, params, iofolder, systime, output)

  # Only print, if output is set to true
  if output == false  return  end

  # Open formatted output file
  open("$iofolder/parameters.dat","w") do f
    # Print header
    println(f,
    "Improved parameters and standard errors for parameterisation of photolysis processes")
    println(f,
    "in the Master Chemical Mechanism (MCM; http://mcm.york.ac.uk/) using:")
    println(f, "\nj / s-1 = l(O3)·(cos(x))^m·exp(-n·sec(x))\n")
    println(f, "with\nl(O3) = l_a0 + l_b0·exp(-O3/l_b1) + l_c0·exp(-O3/l_c1)\n")
    println(f, "Fits and standard errors are derived with Julia package LsqFit v0.6.0")
    println(f, "(https://github.com/JuliaNLSolvers/LsqFit.jl.git).")
    println(f, "created $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))\n")
    println(f,"   l_a0 / s-1           l_b0 / s-1         l_b1 / DU         ",
      "l_c0 / s-1         l_c1 / DU           m              n         Reaction")

    # Loop over reactions
    for i = 1:length(jvals[1].rxn)
      # Print parameters, errors, and reaction label to output file
      @printf(f,"(%6.3f±%.3f)e%d    (%6.3f±%.3f)e%d   %7.2f±%5.2f    (%6.3f±%.3f)e%d   %7.2f±%5.2f    %.3f±%.3f    %.3f±%.3f    %s\n",
      params.l[i][1]/10^jvals[1].order[i], params.sigma[i][1]/10^jvals[1].order[i], Int(jvals[1].order[i]),
      params.l[i][2]/10^jvals[1].order[i], params.sigma[i][2]/10^jvals[1].order[i], Int(jvals[1].order[i]),
      params.l[i][3], params.sigma[i][3],
      params.l[i][4]/10^jvals[1].order[i], params.sigma[i][4]/10^jvals[1].order[i], Int(jvals[1].order[i]),
      params.l[i][5], params.sigma[i][5], params.m[i], params.sigma[i][6], params.n[i], params.sigma[i][7],
      jvals[1].rxn[i])
    end
  end

  # Open csv output file
  open("$iofolder/parameters.csv","w") do f
    # Print header
    println(f, "l_a0,l_b0,l_b1,l_c0,l_c1,m,n,sigma l_a0,sigma l_b0,sigma l_b1,",
      "sigma l_c0,sigma l_c1,sigma m,sigma n,rxn label")

    # Loop over reactions
    for i = 1:length(jvals[1].rxn)
      # Print parameters, errors, and reaction label to output file
      @printf(f,"%.3e,%.3e,%.2f,%.3e,%.2f,%.3f,%.3f,%.3e,%.3e,%.2f,%.3e,%.2f,%.3f,%.3f,%s\n",
      params.l[i][1], params.l[i][2], params.l[i][3], params.l[i][4], params.l[i][5],
      params.m[i], params.n[i], params.sigma[i][1], params.sigma[i][2],
      params.sigma[i][3], params.sigma[i][4], params.sigma[i][5],
      params.sigma[i][6], params.sigma[i][7], jvals[1].rxn[i])
    end
  end
end #function wrt_newparams


"""
    beautify_chem(reactions)

Format `reactions` with LaTeX for nicer titles in plots.
"""
function beautify_chem(reactions::Vector{String})

  # Initialise output
  chem = String[]
  # Loop over reactions and reformat with LaTeX
  for rxn in reactions
    lstr = replace(rxn, "->" => "\\stackrel{h\\nu}{\\longrightarrow}") # format arrows
    lstr = replace(lstr, " + hv" => "") # remove + hv (now above formatted arrows)
    lstr = replace(lstr, r"([A-Za-z)])([0-9]+)" => s"\1_{\2}") # make numbers in chemical formulas subscripts: {\1}_{\2}
    lstr = replace(lstr, r"\(([0-9]+)" => s"(^{\1}") # except for atom's energetic state
    lstr = replace(lstr, "." => "^{.}") # raise radical dots
    lstr = replace(lstr, r"([0-9]+) " => s"\1\\,") # insert halfspace between stoichiometric indices and species
    lstr = replace(lstr, r"=" => "\\!=\\!") # no spaces between double bonds
    lstr = replace(lstr, r"([A-Za-z0-9])-([A-Za-z])" => s"\1\\!-\\!\2") # format dashes
    lstr = replace(lstr, r"C:" => "{\\ddot C}") # format biradical functions
    lstr = replace(lstr, r"CH:" => "{\\ddot C}H") # format biradical functions
    # define and correct exceptions:
    lstr = replace(lstr, " " => "\\ ") # ensure spaces
    # Ensure unslanted font
    lstr = "\\mathrm{"*lstr*"}"

    # Save reformatted reaction to output
    push!(chem, lstr)
  end

  # Return reformatted output
  return latexstring.(chem)
end #function beautify_chem
