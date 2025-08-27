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
#------------------------------------------------------------------------------


# Load libraries and essential RMCSim functions.
library(readxl) #used to import parameter values
#library(pracma)
#library(stringr)
source("RMCSim.R")

unload.model <- function(mName){
  # Construct names of required files and objects from mName.
  dll_name = paste(mName, "_model", sep="")
  dll_file = paste(dll_name, .Platform$dynlib.ext, sep="")
  
  # Unload DLL if it has been loaded.
  dyn.unload(dll_file)
  
}


PBPK_compile <- function(mName = "PBPK_template") {
  # Compile the MCSim PBPK template model.
  compile_model(mName)
}


PBPK_run <- function(mName = "PBPK_template", model.param.filename = NULL, 
                     model.param.sheetname = NULL, exposure.param.filename = NULL, 
                     exposure.param.sheetname = NULL, 
                     data.times = NULL, adj.parms = NULL, 
                     BW.table = NULL, Freef.table = NULL,
                     water.dose.frac = NULL, load.and.unload = TRUE,
                     rtol=1e-8,atol=1e-8, cust_expo = NULL){
  # This function runs a simulation using the compiled MCSim model 
  # "PBPK_template" using model and exposure parameters as described in the 
  # given spreadsheets.
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
  if (load.and.unload == TRUE) load_model(mName)
  
  parms = initParms()
  # Adjust default parameters by importing from given excel file
  all.parms <- load.model.parameters(filename = model.param.filename, sheetname = model.param.sheetname, parms = parms)
  all.exposure.parms <- load.exposure.parameters(filename = exposure.param.filename, sheetname = exposure.param.sheetname, parms = all.parms$model_parms)
  
  chem <- all.parms$chem
  species <- all.parms$species
  sex <- all.parms$sex
  model.pars <- all.parms$model.param
  
  parms = all.exposure.parms$model_parms
  exp.parms = all.exposure.parms$exp.parms
  
  # Adjust any parameters set with new values in adj.parms
  for (qname in names(adj.parms)) {
    if (qname %in% names(parms)) parms[qname] = adj.parms[qname]
    else if (qname %in% names(exp.parms)) exp.parms[qname] = adj.parms[qname]
    else stop(paste0(qname, " is not an adjustable (.model) parameter."), call. = FALSE)
  }
  
  parms <- initParms(parms)
  Y0 = initStates(parms)
  
  # Check for incompatible parameters
  # Only one urinary excretion pathway should be used.
  if ((parms["k_ustc"] != 0.0) & (parms["k_ven_ustc"] != 0.0)){
    stop("Both urinary storage pathways cannot be used at the same time.", call. = FALSE)
  }
  
  # Define times (h) for simulation.
  times = NULL
  if (!is.null(data.times)){
    # list of times provided by user
    times <- sort(unique(c(0, data.times)))
  } else {
    # define list of times based on length of simulation
    sim_hours = all.exposure.parms$sim.days*24
    times = seq(from=0, to=sim_hours, by=0.25)
  }
  
  # Create Forcing functions for changing BW, Free fraction
  Forc = NULL
  
  # Construct BW table
  if (all.exposure.parms$BW_constant == "N"){
    if (is.null(BW.table)) {
      stop("BW is not constant but no BW table provided.", call. = FALSE)
    }
    BW_times = BW.table$times
    BW_out = BW.table$BW
    parms["BW"] = BW_out[1] #assumes that BW at time 0 is the first value of table
    
  } else if (all.exposure.parms$BW_constant == "Y"){
    BW_times = c(0, tail(times,1))
    BW_out = c(parms["BW"], parms["BW"])
    
  }
  
  # Construct free fraction table
  if(all.parms$Free_constant == "N"){
    if (is.null(Freef.table)) {
      stop("Free fraction in plasma is not constant but no table of values provided.", call. = FALSE)
    }
    Freef_times = Freef.table$times
    Freef = Freef.table$Freef
    
  } else if (all.parms$Free_constant == "Y"){
    Freef_times = c(0, tail(times,1))
    Freef = c(parms["F_free"], parms["F_free"])
    
  }
  
  Forc <- list(cbind(times=BW_times, BW_in=BW_out),
               cbind(times=Freef_times, Free_in=Freef))
  
  # If exposure requires dose function, construct the input
  # For example: drinking water or periodic inhalation
  df_dose = NULL
  df_dose_bol_oral = NULL
  df_dose_per_inhal = NULL
  df_dose_step_inhal = NULL
  df_dose_IV = NULL
  df_dose_step_oral = NULL
  if(!is.null(exp.parms)){
    # if the dose is via drinking water, construct the input
    if(all.exposure.parms$water.dose == "Y"){
      # Function to return linear interpolation of the BW at the given time
      BW_fnc <- approxfun(BW_times, BW_out, rule = 2) 
      # rule = 2 means to use value at the closest data extreme if outside the given time interval
      # to use: BW_at_time.t <- BW_fnc(time.t)
      
      # Check that dose fraction vector is valid (and compute if equal doses during day).
      if (all.exposure.parms$water.equal == "equal"){
        # create the water.dose.frac vector with equal doses throughout the day
        water.dose.frac = 1 / exp.parms[["n.doses_water"]]*rep(1, exp.parms[["n.doses_water"]])
      } else if (all.exposure.parms$water.equal == "unequal") {
        # check that the water.dose.frac vector is valid
        if (sum(water.dose.frac) != 1){
          stop("Sum of dose fractions is not equal to 1.", call. = FALSE)
        } else if (any(water.dose.frac > 1)){
          stop("Dose fractions has an element larger than 1.", call. = FALSE)
        } else if (any(water.dose.frac < 0)){
          stop("Dose fractions has an element smaller than 0.", call. = FALSE)
        }
      }
      
      # Daily water ingestion times (24 hour clock)
      t.dose.interval = (exp.parms[["t.final.dose_water"]] - exp.parms[["t.first.dose_water"]])/(exp.parms[["n.doses_water"]]-1)
      tz_daily = seq(from = exp.parms[["t.first.dose_water"]], to = exp.parms[["t.final.dose_water"]], by = t.dose.interval)
      
      # Initialize vector of all water dose times and of all water dose amounts.
      tz = vector("numeric", all.exposure.parms$sim.days * length(tz_daily))
      dz = vector("numeric", all.exposure.parms$sim.days * length(tz_daily))
      
      # Build a vector to contain the simulation times at which water ingestion
      # occurs and a vector that contains the corresponding doses.
      idx = 1
      for (day in seq(1, all.exposure.parms$sim.days, 1)) { # for each day
        # midpoint time each day
        t.midpoint = (day - 1) * 24 + 12 
        
        # mid-day BW
        BW.day = BW_fnc(t.midpoint)
        
        # total dose for the day (mg/d)
        tot.dose.day = BW.day*exp.parms[["dose_water"]] # total dose for the day (mg/d)
        
        # Drinking water dose per event (mg).
        dose.per.event = tot.dose.day*water.dose.frac
        
        #add to vector of times and vector of doses
        tz[seq(idx, idx+length(tz_daily)-1, 1)] = (day - 1) * 24 + tz_daily
        dz[seq(idx, idx+length(tz_daily)-1, 1)] = dose.per.event
        idx = idx + length(tz_daily)
      }
      
      # Data frame containing bolus dosing events.
      df_dose_bol_oral = data.frame(var="A_glumen", 
                           time=tz, value=dz,
                           method="add")
      
      # Data frame containing updates to total input to system.
      df_in = df_dose_bol_oral
      df_in["var"] = "A_in"
      
      # Combine and create sorted data frame.
      df_dose_bol_oral = rbind(df_dose_bol_oral, df_in)
      df_dose_bol_oral = df_dose_bol_oral[order(df_dose_bol_oral$time), ]
      
    } #end if(water.dose == "Y")
  
  if(all.exposure.parms$inhal.dose == "Y"){
    # if the dose is via periodic inhalation, construct the input
    if ("N.days.exp" %in% names(exp.parms)) {
      # Inhaled concentration (ppm).
        dose = parms[["Conc_init"]]
        dose_per_event = c(dose, 0)
        # Daily exposure pattern (time on, time off).
        time.exp.starts = exp.parms[["time.exp.starts"]]
        length.exp.day = exp.parms[["length.exp.day"]]
        N.days.exp = exp.parms[["N.days.exp"]]
        tz_daily = c(time.exp.starts, time.exp.starts+length.exp.day)
      
        # Vector of exposure change times.
        tz = vector("numeric", all.exposure.parms$sim.days * 2)
      
        # Build a vector to contain the simulation times at which exposure changes.
        idx = 1
        for (day in seq(1, all.exposure.parms$sim.days, 1)) {
          tz[idx:(idx+1)] = (day - 1) * 24 + tz_daily
          idx = idx + 2
        }
      
        # Build a vector containing the indexes of simulation times at which exposure changes.
        keep = N.days.exp*length(tz_daily)
        toss = (7-N.days.exp)*length(tz_daily)
        idx = 1
        toss.idx = c()
        for (week in seq(1, all.exposure.parms$sim.days/7, 1)) {
          toss.idx[seq(idx, idx+toss-1, 1)] = seq(keep+1 + (week-1)*7*length(tz_daily), keep+toss + (week-1)*7*length(tz_daily), 1)
          idx = idx + toss
        }
        # Data frame containing inhalation events.
        df_dose_per_inhal = data.frame(var="Conc", time=tz, value=dose_per_event,
                                       method="rep")
        # Remove weekend times.
        df_dose_per_inhal = df_dose_per_inhal[-toss.idx,]
        
        if (!is.na(model.pars$Q_ccinh)) {
          df_q = df_dose_per_inhal
          df_q$var = "Q_cc"
          df_q$value = c(model.pars$Q_ccinh, model.pars$Q_cardiacc)
          df_dose_per_inhal = rbind(df_dose_per_inhal,df_q)
          Y0["Q_cc"] = model.pars$Q_ccinh
          parms["Q_cardiacc"] = model.pars$Q_ccinh
        }
        # Sorted data frame.
        df_dose_per_inhal = df_dose_per_inhal[order(df_dose_per_inhal$time), ]
        
        } else { # end if (!is.null(exp.parms[["N.days.exp"]]))
          # Exposure assumed to start exposure at beginning of simulation
          if (!is.null(exp.parms["inh.stop.time"])) { # Ends before end of simulation
            df_dose_step_inhal = data.frame(var=c("Conc", "Q_cc"),
                                            time=rep(exp.parms["inh.stop.time"],2),
                                            value=c(0, model.pars$Q_cardiacc),
                                            method=rep("rep",2))
            }
          if (!is.na(model.pars$Q_ccinh)) {
            parms["Q_cardiacc"] = model.pars$Q_ccinh
            Y0 = initStates(parms)
          }
        }
    } #end if(inhal.dose == "Y")
  } #end if(!is.null(exp.parms))

  # IV infusion
  if (!is.null(all.exposure.parms$T_iv_infuse)) {
    iv_dose_tot = parms["iv_dose"]*parms["BW"] # iv dose in mg
    # reset initial state of amount in blood to zero
    parms["iv_dose"] = 0 
    parms = initParms(parms)
    Y0 = initStates(parms)
    rate_infused = iv_dose_tot/all.exposure.parms$T_iv_infuse
    Y0["R_IV"] = rate_infused
    df_dose_IV = data.frame(var=c("R_IV"), time=c(all.exposure.parms$T_iv_infuse),
                            value=0, method="rep")
  }
  
  # Continuous Oral Dose
  if(!is.null(all.exposure.parms$R_oral)){
    Y0["R_oral"] = all.exposure.parms$R_oral*parms[["BW"]]/24 # rate of oral dose in mg/h
    if(!is.null(df_dose_per_inhal)) { # If an inhalation dosing schedule was created (also), we'll use that
      df_dose_step_oral = subset.data.frame(df_dose_per_inhal,var=="Conc")
      df_dose_step_oral$var <- "R_oral"
      df_dose_step_oral$value[df_dose_step_oral$value>0] <- Y0["R_oral"]
      } else {
        if (!is.null(all.exposure.parms$T_oral_rate)) {
          df_dose_step_oral = data.frame(var="R_oral", 
                                         time=all.exposure.parms$T_oral_rate,
                                         value=0, method="rep")
          }
      }
  }
    
  df_dose = rbind(df_dose_bol_oral, df_dose_per_inhal, df_dose_step_inhal, 
                  df_dose_IV, df_dose_step_oral, cust_expo)
  
  # Sorted data frame.
  if(!is.null(df_dose)) df_dose = df_dose[order(df_dose$time), ]
  
  # Determine initial states for cases of endogenous production when given by a 
  # venous blood concentration at steady state.
  # Assume no other exposure. **This is taken care of within the 
  # compute_endog_rate function and needs to be updated whenever new parameters 
  # that affect the initial states in the model are added.**
  
  # Compute endogenous production rate and SS values based on initial BW, Free fraction in plasma 
  Forc_SS = NULL
  Forc_SS <- list(cbind(times=c(0, tail(times,1)), 
                        BW_in=c(parms["BW"], parms["BW"])),
               cbind(times=c(0, tail(times,1)), 
                     Free_in=c(parms["F_free"], parms["F_free"])))
  
  if (!is.null(all.exposure.parms$C_ven_SS)) {
    if (all.exposure.parms$C_ven_SS != 0.0) {
      endog_rate = compute_endog_rate(c_data = all.exposure.parms$C_ven_SS, mName = mName, 
                                      parms = parms, Forc = Forc_SS)
      parms["R_0bgli"] = endog_rate
      print(paste0("Calculated Endogenous Production Rate: ", endog_rate))
    }
  }
  # if there is an endogenous rate, 
  # adjust the initial conditions to include what is at steady state
  if (parms["R_0bgli"] != 0.0){ 
    Y0_SS = find_SS_PBPK_noDose(mName = mName, parms = parms, Forc = Forc_SS)
    Y0_SS["V_max_li_t"] = 0.0  # this state should not change due to endog. prod.
    Y0_SS["V_max_om_t"] = 0.0  # this state should not change due to endog. prod.
    Y0 = Y0 + Y0_SS #have original states that reflect dosing plus the states that arise from the background concentration
  }
  
  # Run simulation using the actual exposure information. (Assume units in mg/L)
  out = run_model(mName, times, Y0=Y0, parms=parms, forcing = Forc,
                  event_list=list(data=df_dose),
                  rtol=rtol, atol=atol)
  
  out_df <- as.data.frame(out)
  
  # Only return simulation values at requested time points
  out_df <- out_df[which(out_df$time %in% times),]
  
  out_df$time.days <- out_df$time/24  # Convert hours to days
  out_df$time.hr <- out_df$time
  
  # Unload the dll file unless user specifies not to do so
  if(load.and.unload == TRUE) unload.model(mName)
  
  return(out_df)
  
}

# For use with background concentrations, find the steady state concentrations 
# and the rate of endogenous production
find_SS_PBPK_noDose <- function(mName, parms, Forc){
  
  # Set all dosing parameters to zero so that only endogenous rate is being used
  parms["iv_dose"] = 0.0
  parms["oral_dose_init"] = 0.0
  parms["Conc_init"] = 0.0
  
  # Recalculate parameters and initial states based on new dosing parameters
  parms = initParms(parms)
  Y0 = initStates(parms)
  
  t_data = 1e6*7*24  # initial time value, 1e6 weeks in hours
  delta_SS = Inf
  
  while (delta_SS > 0.001){  # Check if at steady state 
    
    times_SS <- c(0, t_data)
    out_SS1 <- run_model(mName, times_SS, Y0, parms, forcing = Forc)
    out_SS1.df <- as.data.frame(out_SS1)
    
    times_SS <- c(0, t_data*10)
    out_SS2 <- run_model(mName, times_SS, Y0, parms, forcing = Forc)
    out_SS2.df <- as.data.frame(out_SS2)
    
    delta_SS <- abs(out_SS1.df[2,"C_ven"] - out_SS2.df[2,"C_ven"])/out_SS1.df[2,"C_ven"]
    t_data = t_data*1.5
  }
  
  Y0 = out_SS1[nrow(out_SS1), names(Y0)]
  
  return(Y0)
  
}

# Function to compute endogenous rate of production necessary to have 
# concentration c_data in venous blood at steady state
compute_endog_rate <- function(c_data, mName, parms, Forc){
  
  # Set all dosing parameters to zero so that only endogenous rate is being used
  parms["iv_dose"] = 0.0
  parms["oral_dose_init"] = 0.0
  parms["Conc_init"] = 0.0
  
  # Recalculate parameters and initial states based on new dosing parameters
  parms = initParms(parms)
  Y0 = initStates(parms)
  
  delta_SS = Inf
  t_data <- c(1e6*7*24.0) #initial value for when to compute SS
  
  opt.theta.bounds = TRUE # is optimal theta value at search bounds
  lower.bound = 1e-1
  upper.bound = 1e1
  
  while (delta_SS > 0.001){  # check for difference larger than 0.1%
    
    while (opt.theta.bounds){ # if optimal theta value is at search bounds
      
      # Initial guess for the parameter value
      theta_init_all <- c(0.5, 1, 2.5)
      cost_opt = Inf
      
      for(theta_init in theta_init_all){
        
        opt_res = optim(par = theta_init, fn = endog_cost_fun, gr = NULL,
                        t_data = t_data, c_data = c_data, mName = mName, 
                        Y0 = Y0, parms = parms, Forc = Forc, method = "Brent", 
                        lower = lower.bound, upper = upper.bound)  #need to check what to do about the search bounds...
        
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
    # Note, this does not check that the simulation is at the desired 
    # steady state concentration.
    times_SS <- c(0, t_data)
    parms["R_0bgli"] <- theta_opt
    
    out_SS1 <- run_model(mName, times_SS, Y0, parms, forcing = Forc)
    out_SS1 <- as.data.frame(out_SS1)
    
    times_SS <- c(0, t_data*10)
    out_SS2 <- run_model(mName, times_SS, Y0, parms, forcing = Forc)
    out_SS2 <- as.data.frame(out_SS2)
    
    delta_SS <- abs(out_SS1[2,"C_ven"] - out_SS2[2,"C_ven"])/out_SS1[2,"C_ven"]
    t_data = t_data*1.5
  } #end while (delta_SS > 0.01)
  
  return(theta_opt)
}

# Cost function used when computing the endogenous rate of production
endog_cost_fun <- function(theta, mName, t_data, c_data, Y0, parms, Forc, epsilon = 1.0e-12){
  
  # Set new parms values
  parms["R_0bgli"] = theta[1]
  Y0 = initStates(parms)
  
  t_data = c(0, t_data)
  
  # Obtain model concentrations
  out = run_model(mName, t_data, Y0, parms, forcing = Forc)
  c_model = out[2, "C_ven"]
  
  # Compute sum of squared relative errors
  SE = ((c_data - c_model) / (c_data + epsilon))**2
  SSE = sum(SE)
  
  return(SSE)
  
}


# Functions to set parameter values: model and exposure

load.model.parameters <- function(filename, sheetname = NULL, parms){
  #Set parameters based on importing a spreadsheet of information given by filename
  
  # Set location of input spreadsheets and create full filename
  data.loc = "Inputs/"
  completefile <- paste0(data.loc, filename)
  
  # Identify subtables within the spreadsheet based on the keyword "Code"
  det.table.rows <- read_excel(completefile, sheet = sheetname)
  code.rows <- which(det.table.rows$Code == "Code")
  table1.end <- code.rows[1]-1
  table2.start <- code.rows[1]+1
  table2.end <- code.rows[2]-1
  table3.start <- code.rows[2]+1
  table3.end <- code.rows[3]-1
  table4.start <- code.rows[3]+1
  table4.end <- code.rows[4]-1
  table5.start <- code.rows[4]+1
  table5.end <- code.rows[4]+1+3
  
  #Read in first table - model information
  model.info.tib <- read_excel(completefile, range = cell_rows(1:table1.end), sheet = sheetname)
  model.info <- list()
  for (ii in model.info.tib$Code){
    model.info[ii] <- model.info.tib$Value[[which(model.info.tib$Code == ii)]]
  }
  
  #read in second table - model parameter values
  model.param.tib <- read_excel(completefile, range = cell_rows(table2.start:table2.end), sheet = sheetname)
  model.param <- list()
  for (ii in model.param.tib$Code){
    model.param[ii] <- model.param.tib$Value[[which(model.param.tib$Code == ii)]]
  }
  
  #read in third table - compartment parameter values (V, Q, P)
  comp.param.tib <- read_excel(completefile, range = cell_rows(table3.start:table3.end), sheet = sheetname)
  comp.param <- list()
  for (ii in comp.param.tib$Code){
    comp.param[ii] <- comp.param.tib$Value[[which(comp.param.tib$Code == ii)]]
  }
  
  #read in fourth table - other model parameter values (typically not used in .model file)
  other.model.param.tib <- read_excel(completefile, range = cell_rows(table4.start:table4.end), sheet = sheetname)
  other.model.param <- list()
  for (ii in other.model.param.tib$Code){
    other.model.param[ii] <- other.model.param.tib$Value[[which(other.model.param.tib$Code == ii)]]
  }
  
  #read in fifth table - optional plotting parameters
  #Plotting parameters not currently used/implemented!
  plot.info.tib <- read_excel(completefile, range = cell_rows(table5.start:table5.end), sheet = sheetname)
  plot.info <- list()
  for (ii in plot.info.tib$Code){
    plot.info[ii] <- plot.info.tib$Value[[which(plot.info.tib$Code == ii)]]
  }
  
  # Map input params into "parms" list for .model
  # parms = initParms()
  #parms = c()
  
  # !is.null protects from errors due to a parameter not existing in the spreadsheet
  # !is.na checks that the value should be different from the default
  
  # parms not included in MCSim file:
  for (ii in names(parms)){
    if (!is.null(model.param[[ii]])){
      if (!is.na(model.param[[ii]])) parms[ii] <- model.param[[ii]]
    }
    if (!is.null(comp.param[[ii]])){
      if (!is.na(comp.param[[ii]])) parms[ii] <- comp.param[[ii]]
    }
  }
  
  if (!is.null(model.info[["num.blood.comp"]])){
    if (model.info$num.blood.comp == "1") parms["single_blood"] <- 1.0
    if (model.info$num.blood.comp == "2") parms["single_blood"] <- 0.0
  }
  
  if (!is.null(model.info[["venous_ss"]])){
    if (model.info$venous_ss == "Y") parms["venous_ss"] <- 1.0
    if (model.info$venous_ss == "N") parms["venous_ss"] <- 0.0
  }
  
  if (!is.null(model.info[["arterial_ss"]])){
    if (model.info$arterial_ss == "Y") parms["arterial_ss"] <- 1.0
    if (model.info$arterial_ss == "N") parms["arterial_ss"] <- 0.0
  }
  
  if (!is.null(model.info[["exist_lung"]])){
    if (model.info$exist_lung == "Y") parms["exist_lung"] <- 1.0
    if (model.info$exist_lung == "N") parms["exist_lung"] <- 0.0
  }
  
  if (!is.null(model.info[["GE_ss"]])){
    if (model.info$GE_ss == "Y") parms["GE_ss"] <- 1.0
    if (model.info$GE_ss == "N") parms["GE_ss"] <- 0.0
  }
  
  if (!is.null(model.info[["CDSW"]])){
    if (model.info$CDSW == "Y") parms["CDSW"] <- 1.0
    if (model.info$CDSW == "N") parms["CDSW"] <- 0.0
  }
  
  
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
  
  # Constants with units of 1/h (rate constants)
  perT.constants <- c("k_bilec", "k_urinec", "k_fst", "k_fecesc", "k_absgi", 
                      "k_absli", "k_absli2", "k_unabs", "k_ustc", 
                      "k_off_li", "k_met_lic", "k_met_luc", "k_met_omc", "k_enz_resyn")
  for (qname in perT.constants) {
    parms[qname] = parms[qname]/convertT
  }
  
  # Constants with units of L/h (volume flow rate constants)
  VperT.constants <- c("Q_cardiacc", "k_loss")
  for (qname in VperT.constants) {
    parms[qname] = parms[qname]*convertV/convertT
  }
  if ("Q_ccinh" %in% names(model.param)) {
    model.param["Q_ccinh"] = model.param[["Q_ccinh"]]*convertV/convertT
  }
  
  
  # Constants with units of mg/h (maximum metabolic rate, resorption maximum)
  MperT.constants <- c("V_max_reabsc", "V_max_bind_li", "V_max_lic", "V_max_omc",
                       "V_max_luc")
  for (qname in MperT.constants) {
    parms[qname] = parms[qname]*convertM/convertT
  }
  
  # Constants with units of mg/L (affinity constants)
  MperV.constants <- c("K_m_reabs", "K_m_bind_li", "K_m_li", "K_m_om", "K_m_lu")
  for (qname in MperV.constants) {
    parms[qname] = parms[qname]*convertM/convertV
  }
  
  # Constants with units of L/mg (2nd order rate constant for enzyme destruction)
  VperM.constants <- c("k_enz_loss")
  for (qname in VperM.constants) {
    parms[qname] = parms[qname]*convertV/convertM
  }
  
  # C.units gives the units of the concentration data
  # A.units gives the units of the amount data
  # it is assumed that all parameters and dose information is in mg, L
  list_of_parms <- list("model_parms" = parms, sim.days = model.param$sim.time, 
                        C.units = plot.info$C.units, 
                        A.units = plot.info$A.units, species = model.info$species, 
                        sex = model.info$sex, chem = model.info$chem.name, 
                        Free_constant = model.info$Free_constant,
                        model.param = model.param)
  return(list_of_parms)
}

load.exposure.parameters <- function(filename, sheetname = NULL, parms){
  #Set parameters based on importing a spreadsheet of information given by filename
  
  # Set location of input spreadsheets and create full filename
  data.loc = "Inputs/"
  completefile <- paste0(data.loc, filename)
  
  # Identify subtables within the spreadsheet based on the keyword "Code"
  det.table.rows <- read_excel(completefile, sheet = sheetname)
  code.rows <- which(det.table.rows$Code == "Code")
  table1.end <- code.rows[1]-1
  table2.start <- code.rows[1]+1
  table2.end <- code.rows[2]-1
  table3.start <- code.rows[2]+1
  table3.end <- code.rows[2]+1+3
  
  #Read in first table - model information
  model.info.tib <- read_excel(completefile, range = cell_rows(1:table1.end), sheet = sheetname)
  model.info <- list()
  for (ii in model.info.tib$Code){
    model.info[ii] <- model.info.tib$Value[[which(model.info.tib$Code == ii)]]
  }
  
  #Read in second table - exposure parameter values
  model.param.tib <- read_excel(completefile, range = cell_rows(table2.start:table2.end), sheet = sheetname)
  model.param <- list()
  for (ii in model.param.tib$Code){
    model.param[ii] <- model.param.tib$Value[[which(model.param.tib$Code == ii)]]
  }  
  
  # Map input params into "parms" list for .model
  # parms included in MCSim file:
    # !is.na checks that the value should be different from the default
  for (ii in names(parms[names(parms) %in% names(model.param)])){
    if (!is.na(model.param[[ii]])) parms[ii] <- model.param[[ii]]
  }
  
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
  if (model.info$oral.dose == "Y") parms["oral_dose_init"] <- model.param$dose_oral*convertM
  
  T_iv_infuse <- NULL
  if (model.info$IV.dose == "Y") {
    parms["iv_dose"] = model.param$dose_iv*convertM
    # Infusion time for IV dose
    
    if (!is.na(model.param$T_iv_infuse)) T_iv_infuse <- model.param$T_iv_infuse*convertT
  }
  
  # Continuous oral dose rate
  R_oral <- NULL
  if (!is.null(model.param$R_oral)) {
    if (!is.na(model.param$R_oral)) R_oral <- model.param$R_oral
  }
  
  # Time to apply continuous oral dose rate until
  T_oral_rate <- NULL
  if (!is.null(model.param$T_oral_rate)) {
    if (!is.na(model.param$T_oral_rate)) T_oral_rate <- model.param$T_oral_rate*convertT
  }
  
  # Background venous blood concentration
  C_ven_SS = NULL
  if (!is.na(model.param[["C_ven_SS"]])) C_ven_SS <- model.param[["C_ven_SS"]]
  
  # Create structure for parameters for periodic exposures
  exp.parms <- NULL

  if (model.info$inhal.dose == "Y"){
    if (!is.na(model.param$Conc_init)) parms["Conc_init"] <- model.param$Conc_init
    if (!is.na(model.param$NCH)) parms["NCH"] <- model.param$NCH
    if (!is.na(model.param$VCHC)) parms["VCHC"] <- model.param$VCHC*convertV
    if (!is.na(model.param$KL)) parms["KL"] <- model.param$KL/convertT
    if (!is.na(model.param$inh.stop.time)) exp.parms["inh.stop.time"] <- model.param$inh.stop.time*convertT
    if (!is.na(model.param$time.exp.starts)) exp.parms["time.exp.starts"] <- model.param$time.exp.starts*convertT
    if (!is.na(model.param$length.exp.day)) exp.parms["length.exp.day"] <- model.param$length.exp.day*convertT
    if (!is.na(model.param$N.days.exp)) exp.parms["N.days.exp"] <- model.param$N.days.exp
  }
  
  if (model.info$water.dose == "Y"){
    if (!is.na(model.param$n.doses_water)) exp.parms["n.doses_water"] <- model.param$n.doses_water
    if (!is.na(model.param$t.first.dose_water)) exp.parms["t.first.dose_water"] <- model.param$t.first.dose_water
    if (!is.na(model.param["t.final.dose_water"])) exp.parms["t.final.dose_water"] <- model.param["t.final.dose_water"]
    if (!is.na(model.param$t.final.dose_water)) exp.parms["t.final.dose_water"] <- model.param$t.final.dose_water
    if (!is.na(model.param$dose_water)) exp.parms["dose_water"] <- model.param$dose_water
  }
  
  # it is assumed that all parameters and dose information is in mg, L
  list_of_parms <- list("model_parms" = parms, sim.days = model.param$sim.time, 
                        water.dose = model.info$water.dose, 
                        inhal.dose = model.info$inhal.dose,
                        exp.parms = exp.parms,
                        water.equal = model.info$water.equal, 
                        T_iv_infuse = T_iv_infuse, C_ven_SS = C_ven_SS, 
                        BW_constant = model.info$BW_constant, 
                        R_oral = R_oral, T_oral_rate = T_oral_rate)
  return(list_of_parms)
}

# Error analysis functions

perc.diff <- function(model, data, tol=1e-6){
  # Compute percent difference between model and data values relative to data
  perc.diff <- 100*(data-model)*((data-model)>tol)/data
  return(perc.diff)
}

perc.diff.scale <- function(model, data, fig.scale){
  # Compute percent difference between model and data values relative to scale 
  # of the digitized figure
  perc.diff <- 100*((model-data)/fig.scale)
  return(perc.diff)
}

perc.diff.calc <- function(model, data){
  # Compute percent difference between two models 
  # (relative to the average value of the models)
  perc.diff.calc <- 100*((model-data)/((model+data)/2))
  return(perc.diff.calc)
}
