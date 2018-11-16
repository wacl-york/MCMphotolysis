"""
    fit_jold(jvals, sza)

Derive MCM parameterisations for DataFrame `jvals` with `sza`-dependent _j_ values.
    j / s^-1 = l·cos^m(sza)·exp(-n·sec(sza))]

The Vector of solar zenith angles `sza` must be in rad.
"""
function fit_jold(jvals)

  # Initialise arrays for fitting data
  fit = []; sigma = []; rmse = []; R2 = []
  ydata = [jvals.jval[i] ./= 10^jvals.order[i]  for i = 1:length(jvals.jval)]

  # Loop over all j data
  for i = 1:length(jvals.jval) #1:length(jvals)

    # Define fit function with initial guesses
    p0 = [1.35jvals.jval[i][1],0.8,0.3]
    # Fit at order 0 to increase accuracy

    # Derive fit
    push!(fit, curve_fit(jold, jvals.rad, ydata[i], p0))
    # Derive sigma with 95% confidence
    push!(sigma, margin_error(fit[i],0.05))
    # Calculate statistical data for RMSE and R^2
    ss_err = sum(fit[i].resid.^2)
    ss_tot = sum((jvals.jval[i].-mean(jvals.jval[i])).^2)
    # RMSE
    push!(rmse, √(ss_err/fit[i].dof))
    # R^2
    push!(R2, 1. - (ss_err/ss_tot))
  end
  jvals.fit = fit;
  jvals.sigma = sigma;
  jvals.RMSE = rmse;
  jvals.R2 = R2
  jvals.l = Float64[fit[i].param[1] for i = 1:length(fit)]
  jvals.m = Float64[fit[i].param[2] for i = 1:length(fit)]
  jvals.n = Float64[fit[i].param[3] for i = 1:length(fit)]

  ydata = [jvals.jval[i] .*= 10^jvals.order[i]  for i = 1:length(jvals.jval)]

  return jvals
end #function fit_j


"""
    getMCMparams(jvals, o3col)

From dictionary `jvals` with j value related data and vector `o3col`
with ozone column values, return a Matrix with the old l parameters for every
ozone column, and vectors of the m and n parameters for the original MCM photolysis
parameterisations of every photolysis reaction in `jvals`.
"""
function getMCMparams(jvals, o3col)
  fit = []
  iO3 = findlast(o3col .≤ 350)
  j350 = fit_jold(jvals[iO3][:jvals], jvals[iO3][:rad])
  m = [f.param[2] for f in j350[:fit]]
  n = [f.param[3] for f in j350[:fit]]
  l = zeros(Float64, length(m), length(jvals))
  for o3 = 1:length(jvals), jmax = 1:length(m)
    l[jmax, o3] = jvals[o3][:jvals][jmax][1]⋅
                  10^jvals[o3][:order][jmax]/10^jvals[iO3][:order][jmax]
  end

  return l, m, n
  # Loop over initial guesses, otherwise drop low o3col value
end


"""
    fitl(l, o3col)

From vector `l` with vectors of l parameters for every ozone column and every
photolysis reaction  and vector `o3col` with the ozone column values, return
a Matrix with the revised parameters for the ozone column dependent new l parameter.
"""
function fitl(l, o3col, rxn)
  # Initialise
  lpar = []; p0 = [0.,1.,100.,1.,100.]
  o3 = convert.(Float64, o3col)
  # Loop over reactions
  for i = 1:length(l[:,1])
    # Fit l parameter to ozone column dependence
    fit = curve_fit(lnew, o3, l[i,:], p0)
    # Adjust initial guesses for non-convergence
    fit, fail = convergel(o3, l[i,:], fit, "p0", 5,
      ["\033[36mINFO:\033[0m Reaction $i ($(rxn[i])) converged after ",
      " interations."])
    # Drop lowest ozone column values for non-convergence
    fit, fail = convergel(o3, l[i,:], fit, "low", 3,
      ["\033[36mINFO:\033[0m Reaction $i ($(rxn[i])) converged after dropping lowest ",
      " O3 column values."])
    # Drop highest ozone column values for non-convergence
    fit, fail = convergel(o3, l[i,:], fit, "high", 3,
      ["\033[36mINFO:\033[0m Reaction $i ($(rxn[i])) converged after dropping highest ",
        " O3 column values."])
    # Warn, if convergence couldn't be reached
    if fail
      print("\033[95mFitting of parameter l did not converge for reaction ")
      println("$i:\033[0m $(string(rxn[i])).")
    end
    # Save improved l parameters
    push!(lpar, fit.param)
  end

  return lpar
end


"""
    convergel(o3col::Vector{Float64}, lpar::Vector{Float64},
             fit::LsqFit.LsqFitResult, test::String, maxtry::Int64, error_msg::Vector{String})

If `fit` of improved parameter l did not converge, depending on keyword `test`,
try to reach converging by adjusting the initial guesses `"p0"` (`trymax` times),
or dropping up to the `trymax` `"low"`est or `"high"`est ozone column values and
refitting `lpar` (Matrix with l parameters of every reaction for every ozone column)
until convergence is reached and print warning on success.
"""
function convergel(o3col::Vector{Float64}, lpar::Vector{Float64},
         fit::LsqFit.LsqFitResult, test::String, maxtry::Int64, error_msg::Vector{String})
  counter = 0; fail = false
  while !fit.converged
    counter += 1
    if test == "p0"
      fit = curve_fit(lnew, o3col, lpar, fit.param)
    elseif test == "low"
      fit = curve_fit(lnew, o3col[1+counter:end], lpar[1+counter:end,i], fit.param)
    elseif test == "high"
      fit = curve_fit(lnew, o3col[counter:end-counter],
        lpar[counter:end-counter,i], fit.param)
    end
    if fit.converged
      println(error_msg[1],counter,error_msg[2])
    elseif counter == maxtry
      fail = true
      break
    end
  end

  return fit, fail
end


# Define fitting function
"""
    jold(χ,p) / s^-1 = p[1]·cos(χ)^p[2]·exp(-p[3]·sec(χ))]
"""
jold(χ,p) = p[1].*(cos.(χ)).^(p[2]).*exp.(-p[3].⋅sec.(χ))

# Calculate j(χ[, O3col])
"""
    jMCM(χ,O3col,lpar,m,n) / s^-1 = lnew(O3col,lpar)·cos^m(χ)·exp(-n·sec(χ))]
    = (lpar[1] + lpar[2]·exp(-O3/lpar[3]) + lpar[4]·exp(-O3/lpar[5]))·cos^m(χ)·exp(-n·sec(χ))]
"""
jMCM(χ,O3col,lpar,m,n) = lnew(O3col,lpar).⋅
                        (cos.(χ)).^m.⋅exp.(-n.⋅sec.(χ))
# Fit ozone column dependency of the l parameter
"""
    lnew(O3,p) / s^-1 = p[1] + p[2]·exp(-O3/p[3]) + p[4]·exp(-O3/p[5])
"""
lnew(O3,p) = p[1].+p[2].⋅exp.(-O3./p[3]).+p[4].⋅exp.(-O3./p[5])
