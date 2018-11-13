#########################################
###  P U B L I C   F U N C T I O N S  ###
#########################################

"""
    fit_j(jvals, sza)

Derive MCM parameterisations for DataFrame `jvals` with `sza`-dependent _j_ values.
    j / s^-1 = l·cos^m(sza)·exp(-n·sec(sza))]

The Vector of solar zenith angles `sza` must be in rad.
"""
function fit_jold(jvals,sza)

  # Initialise arrays for fitting data
  fit = []; sigma = []; rmse = []; R2 = []

  # Loop over all j data
  @showprogress 1 "fit data ..." for i = 1:length(jvals) #1:length(jvals)

    # Define fit function with initial guesses
    p0 = [jvals[i][1],0.8,0.3]

    # Derive fit
    push!(fit, curve_fit(jold, sza, jvals[i], p0))
    # Derive sigma with 95% confidence
    push!(sigma, margin_error(fit[i],0.05))
    # Calculate statistical data for RMSE and R^2
    ss_err = sum(fit[i].resid.^2)
    ss_tot = sum((jvals[i].-mean(jvals[i])).^2)
    # RMSE
    push!(rmse, √(ss_err/fit[i].dof))
    # R^2
    push!(R2, 1. - (ss_err/ss_tot))
  end

  return Dict(:fit => fit, :σ => sigma, :RMSE => rmse, :R2 => R2)
end #function fit_j


###########################################
###  P R I V A T E   F U N C T I O N S  ###
###########################################

# Define fitting function
jold(x,p) = p[1].*(cos.(x)).^(p[2]).*exp.(-p[3]./cos.(x))

# Calculate j(χ[, O3col])
jMCM(χ,O3col,lpar,m,n) = lnew(O3col,lpar).⋅
                        (cos.(χ)).^m.⋅exp.(-n.⋅sec.(χ))
# Fit ozone column dependency of the l parameter
lnew(x,p) = p[1].+p[2].⋅exp.(-x./p[3]).+p[4].⋅exp.(-x./p[5])
