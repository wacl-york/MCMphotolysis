"""
    plot_jold(jvals::TUVdata, params::PhotData, systime::DateTime, iofolder::String, output::Union{Bool,String})

Plot j values saved in `jvals.jval` and parameterisations derived from parameters
in `params` to file `jvalues.pdf` in the `iofolder` together with the
`systime` of creation, if `output` is set to `true` or `"plot"`.
"""
function plot_jold(jvals::TUVdata, params::PhotData, systime::Dates.DateTime,
  iofolder::String, output::Union{Bool,String})
  # Only print, if output is set to true
  if !(output == true || output == "plot")  return  end

  # Format titles with Latex using header of jvals
  ptitle = set_titles(jvals)

  # Open pdf
  opfile = pdf.PdfPages("$iofolder/jvalues.pdf")

  # Loop over all reactions
  pm.@showprogress 0.1 "plot data..." for i=1:length(params.l)
    # define parameters
    # Calculate parameterised values
    χ = collect(0:π/360:π/2)
    jpar = params.l[i]./10.0^jvals.order[i].⋅cos.(χ).^params.m[i].⋅exp.(-params.n[i].⋅sec.(χ))
    # Load TUV data and parameterisation for plotting
    jplt = pyp.load_PlotData(DataFrame(χ=jvals.deg, j=jvals.jval[i]./10.0^jvals.order[i]),
                            label = "TUV data", pt = "s", lc = "black", lt = "None")
    pplt = pyp.load_PlotData(DataFrame(χ=collect(0:0.5:90), j=jpar), lc = "red", lw = 2,
                            label="MCM Parameterisaton")
    # Plot TUV data

    fig, ax = pyp.plot_data(jplt, pplt, ti = ptitle[i],
    xlabel = "solar zenith angle χ",
    ylabel = "j / 10\$^{$(jvals.order[i])}\$ s\$^{-1}\$",
    xlims=(0,90), ylims=(0,nothing), maj_xticks=15, min_xticks=5, ti_offset=-2,
    legpos="lower left")
    # Plot time stamp
    ax.annotate("created $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))",
          xy=[0.02;0.22], xycoords="axes fraction", ha="left", va="bottom")

    # save plot to temporary png
    opfile.savefig(fig)
    close(fig)
  end

  # close pdf
  opfile.close()
end #function plot_j


"""
    plotj(jvals, params, ptitle, O3col, output, iofolder, systime)

Plot TUV data in `jvals.jval` and parameters from `params` together with the plot
title `ptitle`, the `systime` of creation, and the `O3col` (in DU) to the file
`jvalues.<O3col>.pdf` in the `iofolder`.
"""
function plotj(jvals, params, ptitle, O3col, output, iofolder, systime)

  # Only plot output set in output flag
  if typeof(output) == Bool  return  end
  # Ensure output is an array to be able to loop over it
  if output isa Number  output = [output]  end

  # Loop over ozone columns
  for o3 in output
    # Find index for ozone column in O3col
    iO3 = findfirst(O3col .== o3)
    # Open pdf
    opfile = pdf.PdfPages("$iofolder/jvalues.$(o3)DU.pdf")

    # Loop over all reactions
    pm.@showprogress 0.1 "plot j@$(o3)DU..." for i=1:length(params.l)
      # define parameters
      # Calculate parameterised values
      χ = collect(0:π/360:π/2)
      jpar = jMCM(χ, o3, params.l[i], params.m[i], params.n[i])./10.0^jvals[iO3].order[i]
      # Load TUV data and parameterisation for plotting
      jplt = pyp.load_PlotData(DataFrame(χ=jvals[iO3].deg, j=jvals[iO3].jval[i]./
             10.0^jvals[iO3].order[i]), label = "TUV data @$(o3)DU", pt = "s", lc = "black",
             lt = "None")
      pplt = pyp.load_PlotData(DataFrame(χ=collect(0:0.5:90), j=jpar), lc = "red",
             lw = 2, label="MCM Parameterisaton")
      # Plot TUV data

      fig, ax = pyp.plot_data(jplt, pplt, ti = ptitle[i],
      xlabel = "solar zenith angle χ",
      ylabel = "j / 10\$^{$(jvals[iO3].order[i])}\\,\$s\$^{-1}\$",
      xlims=(0,90), ylims=(0,nothing), maj_xticks=15, min_xticks=5, ti_offset=-2,
      legpos="lower left")
      # Plot time stamp
      ax.annotate("created $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))",
            xy=[0.02;0.22], xycoords="axes fraction", ha="left", va="bottom")

      # save plot to temporary png
      opfile.savefig(fig)
      close(fig)
    end

    # close pdf
    opfile.close()
  end
end #function plotj


"""
    plotl(ldat, params, order, O3col, ptitle, iofolder, systime, output)

If `output` is not `false`, from the TUV calculated l parameters `ldat`, the fit
parameters `params`, the `order` of magnitude for each reaction, the overlying
ozone column array `O3col`, and the beautified reaction labels `ptitle`, generate
plots of the fits versus TUV data for the ozone column dependent l parameter in a
file `iofolder/lpar.pdf` stating the `systime` of creation.
"""
function plotl(ldat, params, order, O3col, ptitle, iofolder, systime, output)

  # Only print, if output is set to true
  if typeof(output) == Bool  return  end

  # Initialise array of plots and define x data
  opfile = pdf.PdfPages("$iofolder/lpar.pdf")
  # Loop over reactions
  pm.@showprogress 0.1 "plot l..." for i = 1:length(ptitle)
    # define parameters
    # Calculate parameterised values
    o3 = collect(O3col[1]:O3col[end])
    lpar = lnew(o3, params.l[i])/10.0^order[i] #lnew is defined in fitTUV.jl
    # Load TUV data and parameterisation for plotting
    ldata  = pyp.load_PlotData(DataFrame(o3=O3col, l=ldat[i,:]./10.0^order[i]),
                 label = "l data", pt = "s", lc = "black", lt = "None")
    lparam = pyp.load_PlotData(DataFrame(o3=o3, l=lpar), lc = "red", lw = 2,
                 label="fit")
    # Plot TUV data

    fig, ax = pyp.plot_data(ldata, lparam, ti = ptitle[i],
      xlabel = "ozone column / DU", ylabel = "l / 10\$^{$(order[i])}\\,\$s\$^{-1}\$",
      xlims=(O3col[1],O3col[end]),ti_offset=-2, legpos="upper right")
    # Plot time stamp
    ax.annotate("created $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))",
          xy=[0.2;0.9], xycoords="axes fraction", ha="left", va="bottom")

    # save plot to temporary png
    opfile.savefig(fig)
    close(fig)
  end

  # compile pngs in single pdf and delete pngs
  opfile.close()
end #function plotl


"""
    write_oldparams(jvals, params, iofolder, systime, output)

Write the parameters in `params` and reaction labels stored in `jvals` to the
formatted file `parameters.dat` and semicolon-separated file `parameters.csv`
in the `iofolder` and state the `systime` of creation, if `output` is not `false`
or `"None"`.
"""
function write_oldparams(jvals, params, iofolder, systime, output) #, stats
  # Only print, if output is set to true
  if output == false || output == "None"  return  end

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
    println(f,"\n                 P a r a m e t e r s")
    println(f,"     l / s-1              m              n         Reaction")

    # Loop over reactions
    for i = 1:length(jvals.rxn)
      # Print parameters, statistical data, and reaction label to output file
      @printf(f,"(%6.3f±%.3f)e%d    %.3f±%.3f    %.3f±%.3f    J(%d): %s\n",
      params.l[i]/10.0^jvals.order[i], params.sigma[i][1]/10.0^jvals.order[i], jvals.order[i],
      params.m[i], params.sigma[i][2], params.n[i], params.sigma[i][3],
      jvals.mcm[i], jvals.rxn[i])
    end
  end

  # Open csv output file
  open("$iofolder/parameters.csv","w") do f
    # Print header
    println(f, "l;m;n;sigma l;sigma m;sigma n;MCM rxn number;TUV rxn number;rxn label")

    # Loop over reactions
    for i = 1:length(jvals.rxn)
      # Print parameters, errors, and reaction label to output file
      @printf(f,"%.3e;%.3f;%.3f;%.3e;%.3f;%.3f;J(%d);%d;%s\n",
      params.l[i], params.m[i], params.n[i], params.sigma[i][1], params.sigma[i][2],
      params.sigma[i][3], jvals.mcm[i], jvals.tuv[i], jvals.rxn[i])
    end
  end
end # function wrt_params


"""
    write_params(jvals, params, iofolder, systime, output)

Write the parameters in `params` and reaction labels stored in `jvals` to the
formatted file `parameters.dat` or semicolon-separated file `parameters.csv`
in the `iofolder`, and state the `systime` of creation, if `output` is not
`false` or `"None"`.
"""
function write_params(jvals, params, iofolder, systime, output)

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
      @printf(f,"(%6.3f±%.3f)e%d    (%6.3f±%.3f)e%d   %7.2f±%5.2f    (%6.3f±%.3f)e%d   %7.2f±%5.2f    %.3f±%.3f    %.3f±%.3f    J(%d): %s\n",
      params.l[i][1]/10.0^jvals[1].order[i], params.sigma[i][1]/10.0^jvals[1].order[i], jvals[1].order[i],
      params.l[i][2]/10.0^jvals[1].order[i], params.sigma[i][2]/10.0^jvals[1].order[i], jvals[1].order[i],
      params.l[i][3], params.sigma[i][3],
      params.l[i][4]/10.0^jvals[1].order[i], params.sigma[i][4]/10.0^jvals[1].order[i], jvals[1].order[i],
      params.l[i][5], params.sigma[i][5], params.m[i], params.sigma[i][6], params.n[i], params.sigma[i][7],
      jvals[1].mcm[i], jvals[1].rxn[i])
    end
  end

  # Open csv output file
  open("$iofolder/parameters.csv","w") do f
    # Print header
    println(f, "l_a0;l_b0;l_b1;l_c0;l_c1;m;n;sigma l_a0;sigma l_b0;sigma l_b1;",
      "sigma l_c0;sigma l_c1;sigma m;sigma n;MCM rxn number;TUV rxn number;rxn label")

    # Loop over reactions
    for i = 1:length(jvals[1].rxn)
      # Print parameters, errors, and reaction label to output file
      @printf(f,"%.3e;%.3e;%.2f;%.3e;%.2f;%.3f;%.3f;%.3e;%.3e;%.2f;%.3e;%.2f;%.3f;%.3f;J(%d);%d;%s\n",
      params.l[i][1], params.l[i][2], params.l[i][3], params.l[i][4], params.l[i][5],
      params.m[i], params.n[i], params.sigma[i][1], params.sigma[i][2],
      params.sigma[i][3], params.sigma[i][4], params.sigma[i][5],
      params.sigma[i][6], params.sigma[i][7], jvals[1].mcm[i], jvals[1].tuv[i], jvals[1].rxn[i])
    end
  end
end #function wrt_newparams


"""
    set_titles(reactions::TUVdata)

Format `reactions.rxn` with LaTeX for nicer titles in plots and add MCM reaction
number from `reactions.mcm`.
"""
function set_titles(reactions::TUVdata)

  # Initialise output
  chem = String[]
  # Loop over reactions and reformat with LaTeX
  for (i, rxn) in enumerate(reactions.rxn)
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
    lstr = "J($(reactions.mcm[i])): \\mathrm{$lstr}"

    # Save reformatted reaction to output
    push!(chem, lstr)
  end

  # Return reformatted output
  return latex.latexstring.(chem)
end #function set_titles
