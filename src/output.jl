"""
    plot_j(jvals, iofolder, fit, systime)

Plot j values saved in dataframe jvals and parameterisations derived from parameters
in fit to file parameters.dat in the iofolder together with the time of creation (time).
"""
function plot_jold(jvals::Dict{Symbol,Any},systime::DateTime,iofolder::String,ifile::String)

  # Format titles with Latex using header of jvals
  ptitle = beautify_chem(names(jvals[:jvals]))

  # Initialise array of plots and define x data
  nj = length(jvals[:fit])
  # jplot = Array{PyCall.PyObject}(nj)
  ofile = "$iofolder/$ifile.pdf"
  opfile = pdf[:PdfPages](ofile)

  # Loop over all reactions
  @showprogress 1 "plot data..." for i=1:nj
    # define parameters
    l = jvals[:fit][i].param[1]
    m = jvals[:fit][i].param[2]
    n = jvals[:fit][i].param[3]
    # Calculate parameterised values
    jpar = l.⋅cos.(jvals[:rad]).^m.⋅exp.(-n.⋅sec.(jvals[:rad]))
    # Load TUV data and parameterisation for plotting
    jplt = pyp.load_PlotData(DataFrame(χ=jvals[:deg], j=jvals[:jvals][i]),
                            label = "TUV data", pt = "s", lc = "black", lt = "None")
    pplt = pyp.load_PlotData(DataFrame(x=jvals[:deg], y=jpar), lc ="red", lw=2,
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
