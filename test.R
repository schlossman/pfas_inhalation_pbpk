# test script to demonstrate ODE solver errors causing GI-related state
# variables to go negative when running steady-state simulations with zero
# exposure (but non-zero endogenous production).

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

# Load functions to run the PBPK model template.
source("run_template_model.R")

model.param.filename = "MeOH_template_parameters_Model.xlsx"
model.param.sheetname = "IRIS_model_rat_2gi"
exposure.param.filename = "MeOH_template_parameters_Exposure.xlsx"
exposure.param.sheetname = "Ward_oral_2500"
data.times = simB5_2500$Time_2gi
rtol=1e-8; atol=1e-8; method="lsoda"

model$updateParms()
# Adjust default parameters by importing from given excel file
mparms <- load.model.parameters(filename=model.param.filename,
                                sheetname=model.param.sheetname, model$parms)
parms <- mparms$parms
chem <- mparms$chem
species <- mparms$species
sex <- mparms$sex
Y0=c()
eparms <- load.exposure.parameters(filename=exposure.param.filename, 
                                   sheetname=exposure.param.sheetname, parms)
parms <- eparms$parms
eoparms <- eparms$other_parms
exp_parms <- eoparms$exp.parms

times = sort(unique(c(0, data.times)))

# Construct BW table
BW_times = c(0, max(times))
BW_out = c(1,1)*parms["BW"]
BW_fnc <- approxfun(BW_times, BW_out, rule = 2)

# Free fraction table
Freef_times = c(0, max(times))
Freef = c(1,1)*parms["F_free"]

# Forcing function for these
Forc <- list(cbind(times=BW_times, BW_in=BW_out),
             cbind(times=Freef_times, Free_in=Freef))

print.noquote(paste("Oral dose is:",parms[["oral_dose_init"]],"mg/kg"))

df_dose=NULL

model$updateParms(parms) # Update model$parms based on parms from above.
model$updateY0()

print.noquote(paste("Venous blood conc from endogenous exposure:",eoparms$C_ven_SS,"mg/L"))

compute_endog_rate(model, c_data=eoparms$C_ven_SS, rtol, atol, method)
print.noquote(paste("Calculated Endogenous Production Rate:", model$parms["R_0bgli"]))

print.noquote(paste("Initial amount in GI lumen:",model$Y0[["A_glumen"]]))

# Then, if there is an endogenous rate (may have been a set input), adjust 
# the initial conditions to include what is at steady state
### find_SS_noDose(model=model, Forc=Forc_SS, rtol=rtol, atol=atol, method=method)

##find_SS_noDose <- function(model, Forc, rtol=1e-8, atol=1e-8, method="lsoda"){
  # For use with background concentrations, find the steady state concentrations 
  # and the rate of endogenous production for 'model' given events df 'Forc'
  
  sparms <- model$parms # Save current set of model parameters
  sY0 <- model$Y0 # Save current initial conditions
  # Set dosing parameters to zero so that only endogenous rate is being used
### update_vals() only changes the values in p with elements of np that have the same name
  p <- update_vals(p=sparms, np=c(iv_dose=0, oral_dose_init=0, Conc_init=0))
  # Recalculate dependent parameters based on new dosing parameters
  model$updateParms(p)
  model$updateY0()
### Check update of oral dose
  print.noquote(paste("Initial amount in GI lumen is now:",model$Y0[["A_glumen"]]))
  t_data = c(0,2400); n=length(t_data)  # initial time value, 100 days in hours
  Forc <- list(cbind(times=c(0, max(t_data)), BW_in=c(1,1)*p["BW"]),
               cbind(times=c(0, max(t_data)), Free_in=c(1,1)*p["F_free"]))
  delta_SS = Inf
  while (delta_SS > 0.001){  # Check if at steady state (< 0.1% change)
    out_SS1 <- as.list(model$runModel(times=t_data, forcings=Forc,
                                      rtol=rtol, atol=atol, method=method)[n,])
    t_data = t_data*2
    out_SS2 <- as.list(model$runModel(times=t_data, forcings=Forc,
                                      rtol=rtol, atol=atol, method=method)[n,])
    delta_SS <- abs(out_SS1$C_ven/out_SS2$C_ven -1)
  }
  
  print.noquote("Simulation resuts (final timepoint) that are negative:")
  print.noquote(out_SS2[out_SS2<0])