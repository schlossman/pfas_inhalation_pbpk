#------------------------------------------------------------------------------
# run_template_model.R
#
# This file contains functions to run the PBPK model template.
#
# Author: Amanda Bernstein, September 2019
#   Updated: February 2021 (note, all changes are now tracked on BitBucket)

# Updates:
# Oct 2020
#  - rename all "PFAS" to "PBPK"
#  - rewrite load.parameters() to be more readable
#  - separate load functions for model parameters and exposure parameters
#
# Feb 2021
#  - remove optimization and sensitivity analysis related code to separate file
#
# Aug-Oct 2025
#  - convert to use MCSimMod instead of older R/MCSim commands
#  - Paul Schlosser, U.S. EPA
#------------------------------------------------------------------------------

# Load libraries and essential MCSimMod and PBPK_run() functions.
library(readxl) # Used to import parameter values
library(MCSimMod) # Load MCSimMod

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

# The following commands creates an object called "template" for use of the
# PBPK Model Template, based on the ODE model encoded in the GNU MCSim model
# language with file name "PBPK_template.model", and a function to run it as
# a PBPK Model Template, using the MCSimMod package by Dustin Kapraun.

template <- createModel("PBPK_template") 

update_vals <- function(p, np, stopifwarned=FALSE){
  # Function expecting two named lists, p (existing parameter/values), np (new 
  # parameter values) and optional 'stopifwarned'. 
  # Elements of p are replaced by elements of np that have corresponding names;
  # i.e., if "name" is the name of both an item in both p and np, then p["name"]
  # is set to np["name"]. The rest of p is unchanged. Applies "as.numeric" to
  # values of np. Warnings are provided for names in np that are not in p but
  # a 'stop' only occurs if stopifwarned=TRUE (default value is FALSE).
  for (n in names(np)){
    if (n%in%names(p)){
      p[n] <- as.numeric(np[n])
    } else { warning(paste0("'",n,"' in list np is not in list p.")) }
  }
  if (stopifwarned & any(!names(np)%in%names(p))) {stop()}
  return(p)
}

PBPK_run <- function(model=template, load=TRUE, 
                     model.param.filename=NULL, model.param.sheetname=NULL, 
                     exposure.param.filename=NULL, exposure.param.sheetname=NULL, 
                     data.times=NULL, adj.parms=NULL, BW.table=NULL,
                     Freef.table=NULL, water.dose.frac=NULL, cust_expo=NULL,
                     rtol=1e-8, atol=1e-8, method="lsoda"){
  # This function runs a simulation using the compiled MCSim model ("model")
  # using model and exposure parameters as described in the input spreadsheets.
  #
  # data.times is a vector of time points (h) at which to return the simulation
  #     output adj.parms is a named vector of .model parameters that will be 
  #     adjusted to have the values given in adj.parms, overwriting the values 
  #     given in the model parameter and exposure parameter spreadsheets.
  # adj.parms is a list of 'par_name = value's that replace those in the spreadsheets.
  # BW.table is a list with two elements: a vector of BWs called BW and a vector
  #     of corresponding times in h called times
  # Freef.table is a list with two elements: a vector of values for the free 
  #     fraction of chemical in plasma called Freef and a vector of corresponding
  #     times in h called times
  # water.dose.frac is a vector listing the fraction of dose ingested at each 
  #     time for bolus water dosing (should sum to 1)
  # load.and.unload is a boolean. When true, the model will automatically be 
  #     loaded at the start of this function and unloaded at the end of this 
  #     function. Set to FALSE to stop automatically loading and unloading the 
  #     model (in which case the model must be loaded prior to running this function).
  # cust_expo is a data frame that allows the user to set up a custom exposure 
  #     pattern. Note, this exposure pattern is applied in addition to any 
  #     exposures defined in the exposure parameter spreadsheet. It is included 
  #     in the events list passed to the ODE solver. It must have the columns 
  #     var, time, value, and method. See the documentation for the R package 
  #     deSolve for more information.
  
  # Load the dll file unless user specifies not to do so
  if (load) model$loadModel()
  
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
  
  # Adjust any parameters set with new values in adj.parms:
  parms <- update_vals(parms, adj.parms, stopifwarned=TRUE)
 
  # Check for incompatible parameters
  # Only one urinary excretion pathway should be used.
  if (parms["k_ustc"]*parms["k_ven_ustc"] != 0){
    stop("Both urinary storage pathways cannot be used at the same time.", call. = FALSE)
  }
  
  # Define times (h) for simulation.
  if (!is.null(data.times)){ times = sort(unique(c(0, data.times))) # list of times provided by user
  } else { times = seq(from=0, to=eoparms$sim.days*24, by=0.25) } # define list of times based on length of simulation
  
  # Create Forcing functions for changing BW, Free fraction
  Forc = NULL
  
  # Construct BW table
  if (eoparms$BW_constant == "N"){
    if (is.null(BW.table)) {
      stop("BW is not constant but no BW table provided.", call. = FALSE)
    }
    BW_times = BW.table$times
    BW_out = BW.table$BW
    parms["BW"] = BW_out[1] #assumes that BW at time 0 is the first value of table
  } else if (eoparms$BW_constant == "Y"){
    BW_times = c(0, max(times))
    BW_out = c(1,1)*parms["BW"]
  }
  BW_fnc <- approxfun(BW_times, BW_out, rule = 2)
  # rule = 2 means to use value at the closest data extreme if outside the given time interval
  # to use: BW_at_time.t <- BW_fnc(time.t)
 
  # Construct free fraction table
  if (mparms$Free_constant){ # load.model.params now converts "Y" to 1, "N" to 0
    Freef_times = c(0, max(times))
    Freef = c(1,1)*parms["F_free"]
  } else { # Assume oparms$Free_constant = 0 ~ "N"
    if (is.null(Freef.table)) {
      stop("Free fraction in plasma is not constant but no table of values provided.", call. = FALSE)
    }
    Freef_times = Freef.table$times
    Freef = Freef.table$Freef
  }
  
  Forc <- list(cbind(times=BW_times, BW_in=BW_out),
               cbind(times=Freef_times, Free_in=Freef))
  
  # Start with null dosing data-frames
  df_dose <- df_dose_bol_oral <- step_oral <- NULL
  df_dose_per_inhal <- df_dose_step_inhal <- df_dose_IV <- NULL

  # If exposure requires dose function, construct the input
  # For example: drinking water or periodic inhalation
  if(!is.null(names(exp_parms))){
    # if the dose is via drinking water, construct the input:
    if(eoparms$water.dose == "Y"){
      # Function to return linear interpolation of the BW at the given time
            nz = exp_parms[["n.doses_water"]]
      # Check that dose fraction vector is valid (and compute if equal doses during day).
      if (eoparms$water.equal == "equal"){
        # create the water.dose.frac vector with equal doses throughout the day
        water.dose.frac = rep(exp_parms[["dose_water"]], ndoses)/ndoses
      } else if (eoparms$water.equal == "unequal") {
        # check that the water.dose.frac vector is valid
        if (sum(water.dose.frac) != 1){
          stop("Sum of dose fractions is not equal to 1.", call. = FALSE)
        } else if (any(water.dose.frac > 1)){
          stop("Dose fractions has an element larger than 1.", call. = FALSE)
        } else if (any(water.dose.frac < 0)){
          stop("Dose fractions has an element smaller than 0.", call. = FALSE)
        }
        water.dose.frac = water.dose.frac*exp_parms[["dose_water"]]
      }
      
      # Daily water ingestion times (24 hour clock)
      tz_daily = seq.int(exp_parms[["t.first.dose_water"]], 
                         exp_parms[["t.final.dose_water"]], length.out=nz)
      # Initialize vector of all water dose times and of all water dose amounts.
      tt=floor((0:(days*nz-1))/nz)*24 # Times at start of each day, repeated
      tz = floor((0:(eoparms$sim.days*nz-1))/nz)*24 + tz_daily
      dz = BW_fnc(tz)*water.dose.frac
      
      # Data frame containing bolus dosing events.
      df_dose_bol_oral = data.frame(var="A_glumen", time=tz, value=dz, method="add")
      
      # Data frame containing updates to total input to system.
      df_in = df_dose_bol_oral
      df_in["var"] = "A_in"
      
      # Combine and create sorted data frame.
      df_dose_bol_oral = rbind(df_dose_bol_oral, df_in)
      df_dose_bol_oral = df_dose_bol_oral[order(df_dose_bol_oral$time), ]
      
    } #end if(eoparms$water.dose == "Y")
  
  if(eoparms$inhal.dose == "Y") {
    # Doses and times for single-day inhalation exposures:
    dose_per_event = c(parms[["Conc_init"]], 0)
    tz = c(0,24*eoparms$sim.days) # Default assumption: continuous exposure from time=0 to end of simulation
    if (!is.null(exp_parms$time.exp.starts)) tz[1] = exp_parms$time.exp.starts
    if (!is.null(exp_parms$inh.stop.time))  tz[2] = tz[1] + exp_parms$inh.stop.time
    if (!is.null(exp_parms$length.exp.day))  tz[2] = tz[1] + exp_parms$length.exp.day
    # if the dose is via periodic inhalation (multi-day exposures), construct the time-sequence:
    if ("N.days.exp" %in% names(exp_parms)) {
      # Build a vector to contain the simulation times at which exposure changes.
      days = 1:eoparms$sim.days-1 # All days starting at zero
      hours <- 24*days[(days%%7)<exp_parms$N.days.exp] # Use remainder function to remove days of week after N.days.exp
      tz = as.vector(sort(c(hours+tz[1],hours+tz[2])))
    }
    # Data frame containing inhalation events (either single- or multi-day)
    df_dose_per_inhal = data.frame(var="Conc", time=tz, value=dose_per_event, method="rep")

    if (!is.null(mparms$Q_ccinh)) { # data-frame to adjust Q_cc if Q_ccinh is set
      df_q = df_dose_per_inhal # Same structure / times as df_dose_per_inhal
      df_q$var = "Q_cc"
      df_q$value = c(mparms$Q_ccinh, parms[["Q_cardiacc"]])
      df_dose_per_inhal = rbind(df_dose_per_inhal,df_q)
      parms["Q_cardiacc"] = mparms$Q_ccinh
      Y0["Q_cc"] = mparms$Q_ccinh
      }
    
    # Sorted data frame.
    df_dose_per_inhal = df_dose_per_inhal[order(df_dose_per_inhal$time), ]
    
    } #end if(eparms$inhal.dose == "Y")
  } #end if(!is.null(exp_parms))

  # IV infusion
  if (!is.null(eoparms$T_iv_infuse)) {
    Y0["R_IV"] = parms[["iv_dose"]]*parms[["BW"]]/eparms$T_iv_infuse # rate in mg/h
    parms["iv_dose"] = 0 # reset initial state of amount in blood to zero
    # Then create data-frame to turn off infusion at specified time
    df_dose_IV = data.frame(var=c("R_IV"), time=eparms$T_iv_infuse, value=0, method="rep")
  }
  
  # Continuous Oral Dose
  if(!is.null(eoparms$R_oral)){
    # Convert rate of oral dose from mg/kg/d to mg/h:
    Y0["R_oral"] = eoparms$R_oral*parms[["BW"]]/24 
    if(!is.null(df_dose_per_inhal)) { 
      # If an inhalation dosing schedule was created (also), assume that oral
      # infusion occurs on the same schedule, use the inhalation frame.
      step_oral = subset.data.frame(df_dose_per_inhal,var=="Conc")
      step_oral$var <- "R_oral"
      ons <- step_oral$value > 0
      step_oral$value[ons] <- eoparms$R_oral*BW_fnc(step_oral$time[ons])/24
      } else { if (!is.null(eoparms$T_oral_rate)) { 
        step_oral = data.frame(var="R_oral", time=eoparms$T_oral_rate, value=0, 
                               method="rep") }
      }
    }
    
  df_dose = rbind(df_dose_bol_oral, df_dose_per_inhal, df_dose_step_inhal, 
                  df_dose_IV, step_oral, cust_expo)
  # Sorted data frame.
  alltimes = times
  if(!is.null(df_dose)) { 
    df_dose = df_dose[order(df_dose$time), ]
    alltimes =sort(unique(c(times,df_dose$time)))
  }
  
  # Determine initial states for cases of endogenous production when given by a 
  # venous blood concentration at steady state.
  # Assume no other exposure. **This is taken care of within the 
  # compute_endog_rate function and needs to be updated whenever new parameters 
  # that affect the initial states in the model are added.**
  
  # Compute endogenous production rate and SS values based on initial BW, Free fraction in plasma 
  model$updateParms(parms)
  model$updateY0(Y0)  
  Forc_SS <- list(cbind(times=c(0, tail(times,1)), BW_in=c(1,1)*parms["BW"]),
                  cbind(times=c(0, tail(times,1)), Free_in=c(1,1)*parms["F_free"]))
  
  if (!is.null(eoparms$C_ven_SS)) {
    if (eoparms$C_ven_SS > 0.0) {
      model$parms["R_0bgli"] = compute_endog_rate(c_data=eoparms$C_ven_SS, 
                                                  model=model, Forc=Forc_SS)
      print(paste("Calculated Endogenous Production Rate:", parms["R_0bgli"]))
    }
  }
  # if there is an endogenous rate, 
  # adjust the initial conditions to include what is at steady state
  if (model$parms["R_0bgli"] > 0.0){ 
    Y0_SS = find_SS_noDose(model=model, Forc=Forc_SS)
    model$updateY0(Y0_SS[Y0_SS>0])
      # Revise initial states with the states that arise from the background
      # exposure that are > 0.
  }
  
  # Run simulation using the actual exposure information. (Assume units in mg/L)
  out = model$runModel(times=alltimes, forcings = Forc, events=list(data=df_dose),
                  rtol=rtol, atol=atol, method=method)
    # Only return simulation values at requested time points
  out_df <- as.data.frame(out[out[,"time"]%in%times,])
  out_df[,"time.days"] <- out_df[,"time"]/24  # Convert hours to days
  return(out_df)
}

find_SS_noDose <- function(model, Forc){
# For use with background concentrations, find the steady state concentrations 
# and the rate of endogenous production for 'model' given events df 'Forc'
  
  sparms <- model$parms # Save current set of model parameters
  # Set dosing parameters to zero so that only endogenous rate is being used
  p <- update_vals(p=sparms, np=c(iv_dose=0, oral_dose_init=0, Conc_init=0))
  # Recalculate dependent parameters based on new dosing parameters
  model$updateParms(p)
  
  t_data = c(0, 100*7*24)  # initial time value, 100 weeks in hours
  delta_SS = Inf
  while (delta_SS > 0.001){  # Check if at steady state (< 0.1% change)
    out_SS1 <- model$runModel(times=t_data, forcings=Forc)
    out_SS2 <- model$runModel(times=2*t_data, forcings=Forc)
    delta_SS <- abs(out_SS1[2,"C_ven"]/out_SS2[2,"C_ven"] -1)
    t_data = t_data*2
  }
  model$updateParms(sparms) # Restore model parameters using saved values
  return(out_SS1[-1, names(model$Y0)])
}

compute_endog_rate <- function(model, c_data, Forc){
# Function to compute endogenous rate of production necessary in 'model' to have 
# concentration 'c_data' in venous blood at steady state given events df 'Forc'
  
  sparms <- model$parms # Save current set of model parameters
  # Set dosing parameters to zero so that only endogenous rate is being used
  p <- update_vals(p=sparms, np=c(iv_dose=0, oral_dose_init=0, Conc_init=0))
  
  # Recalculate dependent parameters based on new dosing parameters
  model$updateParms(p)
  model$updateY0()

  delta_SS = Inf
  t_data <- c(0, 100*7*24) #initial value for time to reach SS, 100 weeks in h
  
  opt.theta.bounds = TRUE # is optimal theta value at search bounds
  lower.bound = 1e-1
  upper.bound = 1e1
  
  while (delta_SS > 0.001){  # check for difference larger than 0.1%
    while (opt.theta.bounds){ # if optimal theta value is at search bounds
      # Initial guess for the parameter value
      theta_init_all <- c(0.5, 1, 2.5)
      cost_opt = Inf
      
      for(theta_init in theta_init_all){
        opt_res = optim(par=theta_init, fn=endog_cost_fun, gr=NULL,
                        t_data=t_data, c_data=c_data, model=model, Forc=Forc,
                        method="Brent", lower=lower.bound, upper=upper.bound)  
            # Need to check what to do about the search bounds...
        
        #print("Cost for optimal parameter is ")
        #print(opt_res$value)
        #print("Optimal parameter is ")
        #print(opt_res$par)
        
        if (opt_res$value < cost_opt) {
          cost_opt = opt_res$value
          theta_opt = opt_res$par
        } #end if(opt_res$value < cost_opt)
      } #end for(theta_init in theta_init_all)
      
      if (abs(theta_opt-lower.bound) < 1e-3){
        #print("theta is lower bound")
        lower.bound = lower.bound*0.1
      } else if (abs(theta_opt-upper.bound) < 1e-3){
        #print("theta is upper bound")
        upper.bound = upper.bound*10
      } else {
        #print("theta is neither bound")
        opt.theta.bounds = FALSE
      }
      
    } #end while (opt.theta.bounds)
    
    # Check if at steady state 
    # Note, this does not check that the simulation is at the *target* 
    # steady state concentration, only that it's at steady state.
    p["R_0bgli"] <- theta_opt
    model$updateParms(p)
    model$updateY0()
    
    out_SS1 <- as.data.frame(model$runModel(t_data, forcings=Forc))
    t_data <- t_data*2
    out_SS2 <- as.data.frame(model$runModel(t_data, forcings = Forc))
    
    delta_SS <- abs(out_SS1[2,"C_ven"] - out_SS2[2,"C_ven"])/out_SS1[2,"C_ven"]
    t_data = t_data*2
  } #end while (delta_SS > 0.01)
  
  return(theta_opt)
}

endog_cost_fun <- function(theta, model, t_data, c_data, Forc, epsilon=1.0e-12){
# Cost function used when computing the endogenous rate of production
  
  model$parms["R_0bgli"] = theta[1]  # Set new parms values
  model$updateY0() # Update initial conditions
  out = model$runModel(t_data, forcings=Forc) # Obtain model concentrations

  # Compute sum of squared relative errors
  SSE = sum(((c_data - out[-1, "C_ven"]) / (c_data + epsilon))**2)
  return(SSE)
}


load.model.parameters <- function(filename, sheetname = NULL, parms){
# Functions to set parameter values: model and exposure
  ptables <- read_excel(paste0("Inputs/", filename), sheet = sheetname)
  notNA <- which(!is.na(ptables$Value)) # Identify rows with Value != NA
  ptable <- as.list(ptables$Value[notNA])
  names(ptable) <- ptables$Code[notNA]

  # Identify subtables within the spreadsheet based on the keyword "Code"
  codeNs <- which(names(ptable)=="Code")
 
  model.info <- ptable[1:(codeNs[1]-1)] # Model information values
  
  model.param <- ptable[c((codeNs[1]+1):(codeNs[2]-1),(codeNs[2]+1):(codeNs[3]-1))] # Model parameters ...
    # Chemical-specific parameters are codeNs[1]+1 to codeNs[2]-1 and
    # compartment parameters (V, Q, P) are codeNs[2]+1 to codeNs[3]
  model.param[]=as.numeric(model.param) # Convert values from text to numeric
  parms <- update_vals(parms, model.param[names(model.param)!="Q_ccinh"])
  
  # Other parameter values (typically not used in .model file):
  other.model.param <- ptable[(codeNs[3]+1):(codeNs[4]-1)]
  
  if (!is.null(model.info[["num.blood.comp"]])){
    # Set parms["single_blood"] = 0 if num.blood.comp = "2", otherwise 1
    if (model.info[["num.blood.comp"]]=="2") { 
      parms["single_blood"] = 0 } else parms["single_blood"] = 1
  }
  
  model.info[w1 <- model.info=="Y"] = 1 # Change "Y" & "1" values to 1, remember which those are
  model.info[w0 <- model.info=="N"] = 0 # Change "N" values to 0, remember which those are
  parms <- update_vals(parms, model.info[(w1|w0)&(names(model.info)!="Free_constant")]) # Update model parms with these

  # Unit conversions
  # Determine unit change parameters (model is assumed to run in mg, L, h)
  # Note, BW is still assumed to be in kg
  # convert is multiplied by the provided units to return [mg, L, h]
  convertM <- 1 # Data units in mg
  convertV <- 1 # Data units in L
  convertT <- 1 # Data units in h
  if (model.info$M.units == "ug") convertM <- 1/1e3
  if (model.info$M.units == "ng") convertM <- 1/1e6
  if (model.info$V.units == "mL") convertV <- 1/1e3
  if (model.info$T.units == "min") convertT <- 1/60
  if (model.info$T.units == "days") convertT <- 24
  
  # Constants with units of 1/h (rate constants) are divided by convertT
  parms <- update_vals(parms, 
                       parms[c("k_bilec", "k_urinec", "k_fst", "k_fecesc", 
                               "k_absgi", "k_absli", "k_absli2", "k_unabs", 
                               "k_ustc", "k_off_li", "k_met_lic", "k_met_luc",  
                               "k_met_omc", "k_enz_resyn")]/convertT) 
  
  # Constants with units of L/h (volume flow rate constants) multiplied by convertV/convertT:
  parms <- update_vals(parms, parms[c("Q_cardiacc", "k_loss")]*convertV/convertT) 
  
  # Q_ccinh is used to adjust Q_cc via an input dataframe. It was introduced
  # for the PFAS inhalation modeling and may not be in all input spreadsheets.
  Q_ccinh = NULL
  if ("Q_ccinh" %in% names(model.param)) { 
    Q_ccinh = model.param[["Q_ccinh"]]*convertV/convertT
  }
  
  # Constants with units of mg/h (maximum metabolic rate, resorption maximum) multiplied by convertM/convertT:
  parms <- update_vals(parms, parms[c("V_max_reabsc", "V_max_bind_li", "V_max_lic", 
                                      "V_max_omc","V_max_luc")]*convertM/convertT)
  
  # Constants with units of mg/L (affinity constants) multiplied by convertM/convertV:
  parms <- update_vals(parms, parms[c("K_m_reabs", "K_m_bind_li", "K_m_li", 
                                      "K_m_om", "K_m_lu")]*convertM/convertV)
  
  # k_enz_loss (L/mg) (2nd order rate constant for enzyme destruction) multiplied by convertV/convertM
  parms <- update_vals(parms, parms["k_enz_loss"]*convertV/convertM)
  
  # C.units gives the units of the concentration data
  # A.units gives the units of the amount data
  # it is assumed that all parameters and dose information is in mg, L
  return(list(parms=parms,
              C.units=paste0(model.info$M.units,"/",model.info$V.units),
              A.units=model.info$M.units, species=model.info$species, 
              sex=model.info$sex, chem=model.info$chem.name,
              Q_ccinh=Q_ccinh, Free_constant=model.info$Free_constant))
}

load.exposure.parameters <- function(filename, sheetname = NULL, parms){
  #Set parameters based on importing a spreadsheet of information given by filename
  ptables <- read_excel(paste0("Inputs/", filename), sheet=sheetname)
  notNA <- which(!is.na(ptables$Value)) # Identify rows with Value != NA
  ptable <- as.list(ptables$Value[notNA])
  names(ptable) <- ptables$Code[notNA]
  
  # Identify subtables within the spreadsheet based on the keyword "Code"
  codeNs <- which(names(ptable)=="Code")
  if (length(codeNs)==1) codeNs <- c(codeNs,length(ptable)+1) # If only one "Code" separator, extend to end of table
  
  model.info <- ptable[1:(codeNs[1]-1)] # Exposure information values

  model.param <- ptable[(codeNs[1]+1):(codeNs[2]-1)] # Exposure parameters ...
  model.param[]=as.numeric(model.param) # Convert values from text to numeric
  parms <- update_vals(parms,model.param[names(model.param)%in%names(parms)])
  
  # Unit conversions
  # Determine unit change parameters (model is assumed to run in mg, L, h)
  # Note, BW is still assumed to be in kg
  # convert is multiplied by the provided units to return [mg, L, h]
  convertM <- 1 # Data units in mg
  convertV <- 1 # Data units in L
  convertT <- 1 # Data units in h
  if (model.info$M.units == "ug") convertM <- 1/1e3
  if (model.info$M.units == "ng") convertM <- 1/1e6
  if (model.info$V.units == "mL") convertV <- 1/1e3
  if (model.info$T.units == "min") convertT <- 1/60
  if (model.info$T.units == "days") convertT <- 24
  
  # Bolus oral dosing information
  if (model.info$oral.dose == "Y") parms["oral_dose_init"] = model.param$dose_oral*convertM
  
  T_iv_infuse <- NULL
  if (model.info$IV.dose == "Y") { 
    parms["iv_dose"] = model.param$dose_iv*convertM
    if ("T_iv_infuse"%in%names(model.param)) T_iv_infuse <- model.param$T_iv_infuse*convertT
  }
  
  # Continuous oral dose rate
  R_oral <- NULL
  if ("R_oral"%in%names(model.param)) R_oral <- model.param$R_oral*convertM
  
  # Time to apply continuous oral dose rate until
  T_oral_rate <- NULL
  if ("T_oral_rate"%in%names(model.param)) T_oral_rate <- model.param$T_oral_rate*convertT

  # Background venous blood concentration
  C_ven_SS = NULL
  if ("C_ven_SS"%in%names(model.param)) C_ven_SS <- model.param$C_ven_SS
  
  # Create structure for parameters for periodic exposures
  exp.parms <- list()

  if (model.info$inhal.dose == "Y"){
    if ("Conc_init"%in%names(model.param)) parms["Conc_init"]=model.param$Conc_init
    if ("NCH"%in%names(model.param)) parms["NCH"]=model.param$NCH
    if ("VCHC"%in%names(model.param)) parms["VCH"]=model.param$VCHC*convertV
    if ("KL"%in%names(model.param)) parms["KL"]=model.param$KL/convertT
    if ("inh.stop.time"%in%names(model.param)) exp.parms$inh.stop.time <- model.param$inh.stop.time*convertT
    if ("time.exp.starts"%in%names(model.param)) exp.parms$time.exp.starts <- model.param$time.exp.starts*convertT
    if ("length.exp.day"%in%names(model.param)) exp.parms$length.exp.day <- model.param$length.exp.day*convertT
    if ("N.days.exp"%in%names(model.param)) exp.parms$N.days.exp <- model.param$N.days.exp
  }
  
  if (model.info$water.dose == "Y"){
    if ("n.doses_water"%in%names(model.param)) exp.parms$n.doses_water <- model.param$n.doses_water
    if ("t.first.dose_water"%in%names(model.param)) exp.parms$t.first.dose_water <- model.param$t.first.dose_water
    if ("t.final.dose_water"%in%names(model.param)) exp.parms$t.final.dose_water <- model.param$t.final.dose_water
    if ("dose_water"%in%names(model.param)) exp.parms$dose_water <- model.param$dose_water
  }
  
  # it is assumed that all parameters and dose information is in mg, L
  other_parms <- list(sim.days = as.numeric(model.param$sim.time), 
                        R_oral=R_oral, water.dose = model.info$water.dose, 
                        inhal.dose = model.info$inhal.dose,
                        exp.parms = exp.parms, T_oral_rate = T_oral_rate,
                        water.equal = model.info$water.equal, 
                        T_iv_infuse = T_iv_infuse, C_ven_SS = C_ven_SS, 
                        BW_constant = model.info$BW_constant)
  return(list(parms=parms,other_parms=other_parms))
}

# Error analysis functions
perc.diff <- function(model, data, tol=1e-6){
  # Compute percent difference between model and data values relative to data
  perc.diff <- 100*(data-model)*((data-model)>tol)/data
  return(perc.diff)
 }
# 
perc.diff.scale <- function(model, data, fig.scale){
  # Compute percent difference between model and data values relative to scale 
  # of the digitized figure
  perc.diff <- 100*((model-data)/fig.scale)
  return(perc.diff)
}
# 
perc.diff.calc <- function(model, data){
  # Compute percent difference between two models 
  # (relative to the average value of the models)
  perc.diff.calc <- 100*((model-data)/((model+data)/2))
  return(perc.diff.calc)
}
