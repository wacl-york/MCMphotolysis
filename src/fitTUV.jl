"""
    fit_jold(jvals, sza)

Derive MCM parameterisations for DataFrame `jvals` with `sza`-dependent _j_ values.
    j / s^-1 = l·cos^m(sza)·exp(-n·sec(sza))]

The Vector of solar zenith angles `sza` must be in rad.
"""
function fit_jold(jvals)

  # Initialise arrays for fitting data
  fit = []; sigma = []; rmse = []; R2 = []

  # Loop over all j data
  for i = 1:length(jvals.jval) #1:length(jvals)

    # Define fit function with initial guesses
    p0 = [1.35/10^jvals.order[i].*jvals.jval[i][1],0.8,0.3]
    # Fit at order 0 to increase accuracy

    # Derive fit
    push!(fit, curve_fit(jold, jvals.rad, jvals.jval[i] ./ 10^jvals.order[i], p0))
    # Derive sigma with 95% confidence
    push!(sigma, standard_error(fit[i]))
    # Calculate statistical data for RMSE and R^2
    ss_err = sum(fit[i].resid.^2)
    ss_tot = sum((jvals.jval[i].-mean(jvals.jval[i])).^2)
    # RMSE
    push!(rmse, √(ss_err/fit[i].dof))
    # R^2
    push!(R2, 1. - (ss_err/ss_tot))
  end
  conv = [f.converged for f in fit]
  l = [fit[i].param[1].*10^jvals.order[i] for i = 1:length(fit)]
  m = [fit[i].param[2] for i = 1:length(fit)]
  n = [fit[i].param[3] for i = 1:length(fit)]
  sigma = [10^jvals.order[i].*sigma[i] for i = 1:length(sigma)]
  rmse  = [10^jvals.order[i].*rmse[i] for i = 1:length(rmse)]

  return PhotData(l, m, n, sigma, conv), StatData(rmse, R2)
end #function fit_j


"""
    getMCMparams(jvals, o3col)

From dictionary `jvals` with j value related data and vector `o3col`
with ozone column values, return a Matrix with the old l parameters for every
ozone column, and vectors of the m and n parameters for the original MCM photolysis
parameterisations of every photolysis reaction in `jvals`.
"""
function getMCMparams(jvals, O3col)
  fit = []
  iO3 = findlast(O3col .≤ 350)
  params350, stats350 = fit_jold(jvals[iO3])
  l = zeros(Float64, length(jvals[iO3].rxn), length(O3col))
  for o3 = 1:length(O3col), jmax = 1:length(jvals[iO3].rxn)
    l[jmax, o3] = jvals[o3].jval[jmax][1]/exp(-params350.n[jmax])
  end

  return l, params350
  # Loop over initial guesses, otherwise drop low o3col value
end


"""
    fitl(l, o3col)

From vector `l` with vectors of l parameters for every ozone column and every
photolysis reaction  and vector `o3col` with the ozone column values, return
a Matrix with the revised parameters for the ozone column dependent new l parameter.
"""
function fitl(ldata, order, o3col, rxn)
  # Initialise
  lpar = []; fits = []
  o3 = convert.(Float64, o3col)
  ldata ./= [10^o for o in order]
  # Loop over reactions
  for i = 1:length(ldata[:,1])
    p0 = [0.,ldata[i,1],100.,ldata[i,1],100.]
    # Fit l parameter to ozone column dependence
    fit = curve_fit(lnew, o3, ldata[i,:], p0)
    # Adjust initial guesses for non-convergence
    fit, fits, fail = convergel(o3, ldata[i,:], fit, fits, "p0", 5,
      ["\033[36mINFO:\033[0m Reaction $i ($(rxn[i])) converged after ",
      " interations."])
    # Drop lowest ozone column values for non-convergence
    fit, fits, fail = convergel(o3, ldata[i,:], fit, fits, "low", 3,
      ["\033[36mINFO:\033[0m Reaction $i ($(rxn[i])) converged after dropping lowest ",
      " O3 column values."])
    # Drop highest ozone column values for non-convergence
    fit, fits, fail = convergel(o3, ldata[i,:], fit, fits, "high", 3,
      ["\033[36mINFO:\033[0m Reaction $i ($(rxn[i])) converged after dropping highest ",
        " O3 column values."])
    # Warn, if convergence couldn't be reached
    if fail
      print("\033[95mFitting of parameter l did not converge for reaction ")
      println("$i:\033[0m $(string(rxn[i])).")
    end

    # Calculate errors

    # Save improved l parameters
    l  = deepcopy(fit.param)
    l[[1,2,4]] *= 10^order[i]
    push!(lpar, l)
    push!(fits, fit)
  end
  ldata .*= [10^o for o in order]

  return lpar, fits
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
         fit::LsqFit.LsqFitResult, fits, test::String, maxtry::Int64, error_msg::Vector{String})
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

  return fit, fits, fail
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
