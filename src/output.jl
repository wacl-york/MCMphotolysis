"""
    plot_j(jvals::Dict{Symbol,Any},systime::DateTime,iofolder::String,ifile::String)

Plot j values saved in `jvals[:jvals]` and parameterisations derived from parameters
in `jvals[:fit]` to file `ifile.pdf` in the `iofolder` together with the
`systime` of creation.
"""
function plot_jold(jvals::Dict{Symbol,Any},systime::DateTime,iofolder::String,ifile::String)

  # Format titles with Latex using header of jvals
  ptitle = beautify_chem(names(jvals[:jvals]))

  # Initialise array of plots and define x data
  nj = length(jvals[:fit])
  ofile = "$iofolder/$ifile.pdf"
  opfile = pdf[:PdfPages](ofile)

  # Loop over all reactions
  @showprogress 1 "plot data..." for i=1:nj
    # define parameters
    l = jvals[:fit][i].param[1]
    m = jvals[:fit][i].param[2]
    n = jvals[:fit][i].param[3]
    # Calculate parameterised values
    χ = collect(0:π/360:π/2)
    jpar = l.⋅cos.(χ).^m.⋅exp.(-n.⋅sec.(χ))
    # Load TUV data and parameterisation for plotting
    jplt = pyp.load_PlotData(DataFrame(χ=jvals[:deg], j=jvals[:jvals][i]),
                            label = "TUV data", pt = "s", lc = "black", lt = "None")
    pplt = pyp.load_PlotData(DataFrame(χ=collect(0:0.5:90), j=jpar), lc ="red", lw=2,
                            label="MCM Parameterisaton")
    # Plot TUV data
    fig, ax = pyp.plot_data(jplt, pplt, ti = ptitle[i],
              xlabel = "solar zenith angle χ",
              ylabel = "j / 10\$^{$(Int(jvals[:order][i]))}\$ s\$^{-1}\$",
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
    wrt_params(jvals, iofolder, systime)

Write the parameters, statistical data, and reaction labels stored in the
dictionary `jvals` to the file `parameters.dat` in the `iofolder` and state
the `systime` of creation.
"""
function wrt_params(jvals, iofolder, systime)

  #transform dataframe column symbols to strings
  rxn = string.(names(jvals[:jvals]))

  # Open output file
  open("$iofolder/parameters.dat","w") do f
    # Print header
    println(f,
    "Parameters and statistical data for parameterisation of photolysis processes")
    println(f,
    "in the Master Chemical Mechanism (MCM; http://mcm.york.ac.uk/) using:")
    println(f, "\nj / s-1 = l·(cos(x))^m·exp(-n·sec(x))")
    println(f, "\n\ncreated $(Dates.format(systime,"dd.mm.yyyy, HH:MM:SS"))")
    println(f,"\n                 P a r a m e t e r s               S t a t i s t i c s")
    println(f,"     l / s-1              m              n         RMSE / s-1    R^2      Reaction")

    # Loop over reactions
    for i = 1:length(jvals[:fit])
      # Print parameters, statistical data, and reaction label to output file
      @printf(f,"(%6.3f±%.3f)e%d    %.3f±%.3f    %.3f±%.3f    %.3e    %.4f    %s\n",
      jvals[:fit][i].param[1], jvals[:σ][i][1], Int(jvals[:order][i]),
      jvals[:fit][i].param[2], jvals[:σ][i][2], jvals[:fit][i].param[3], jvals[:σ][i][3],
      jvals[:RMSE][i]⋅10^jvals[:order][i], jvals[:R2][i], rxn[i])
    end
  end
end # function wrt_params


"""
    beautify_chem(reactions)

Format `reactions` with LaTeX for nicer titles in plots.
"""
function beautify_chem(reactions::Vector{Symbol})

  # Initialise output
  chem = String[]
  # Loop over reactions and reformat with LaTeX
  for rxn in strip.(string.(reactions))
    println(rxn)
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
