# model scripts to perform model parameter optimizations, sensitivity analysis, 
# etc.

# Initial versions of functions copied from the run_template_model.R file
#  Created: February 2021

calc_sens <- function(rel_step = 1.0e-8, mName = "PBPK_templateV1", param.filename, sheetname){
  #calculate sensitivity indices for all parameters listed in q
  
  # Load the dll file
  load_model(mName)
  
  # Adjust default parameters
  #data.loc <- data.location(NULL)
  #all.parms = import_params(filename = param.filename, data.loc = data.loc, sheetname = sheetname)
  all.parms <- load.parameters(mName = "PBPK_templateV1", param.filename = param.filename, sheetname = sheetname)
  q <- all.parms$model_parms
  
  #run simulation using parameter values in q
  out_q <- run_sim(mName = "PBPK_templateV1", parms = q, sim.days = all.parms$sim.days)
  
  # Compute AUC for the blood concentration curve
  AUC_q <- trapz(out_q$time, out_q$C_bl)
  # Compute peak blood concentration
  peak_q <- max(out_q[ , "C_bl"])
  
  # Create a dataframe to store sensitivity indices.
  len_q = length(q)
  sens_df = data.frame(param_val=q,
                       S_AUC=double(len_q),
                       Sn_AUC=double(len_q),
                       S_peak=double(len_q),
                       Sn_peak=double(len_q))
  
  # For each parameter in q, compute a local sensitivity indices for both AUC
  # and peak blood concentration.
  for (qname in names(q)) {
    # Perturb this parameter by the "rel_step" amount.
    qx = q
    #qx[qname] = q[qname] * (1 + rel_step)
    qx[qname] = max(rel_step, q[qname] * (1 + rel_step))
    
    # Run a simulation
    out_qx = run_sim(mName = "PBPK_templateV1", parms = qx, sim.days = all.parms$sim.days)
    
    # Compute AUC for the blood concentration curve
    AUC_qx = trapz(out_qx$time, out_qx$C_bl)
    
    # Compute peak blood concentration
    peak_qx = max(out_qx[ , "C_bl"])
    
    # Compute sensitivity indices for AUC.
    sens_df[qname, "S_AUC"] = (AUC_qx - AUC_q) / (qx[qname] - q[qname])
    sens_df[qname, "Sn_AUC"] = sens_df[qname, "S_AUC"] * (q[qname] / AUC_q)
    
    # Compute sensitivity indices for peak.
    sens_df[qname, "S_peak"] = (peak_qx - peak_q) / (qx[qname] - q[qname])
    sens_df[qname, "Sn_peak"] = sens_df[qname, "S_peak"] * (q[qname] / peak_q)
  }
  
  # Return the dataframe containing the sensitivity analysis results.
  return(sens_df)
  
}

sens1 <- function(param, mName = "PBPK_template", param.filename, sheetname, deltaP = 1e-4){
  #NOTE: deltaP is the percent change in param, so the param increases by deltaP*p
  # Performs 1 at a time sensitivity analysis on a given parameter
  
  # Load the dll file
  load_model(mName)
  
  # Adjust default parameters
  data.loc <- data.location(NULL)
  all.parms = import_params(filename = param.filename, data.loc = data.loc, sheetname = sheetname)
  
  # Define times for simulation.
  sim_hours = all.parms$sim.days*24
  times = seq(from=0, to=sim_hours, by=0.1)
  
  Y0 = initStates(all.parms$model_parms)
  
  # Run simulation. (Assume original units in mg/L)
  out0 = run_model(mName, times, Y0, all.parms$model_parms)
  out_df0 <- as.data.frame(out0)
  out_df0$time <- out_df0$time/24  # Convert hours to days
  
  
  # Increase the given parameter by 1%
  param0 <- all.parms$model_parms[param]
  all.parms$model_parms[param] = all.parms$model_parms[param]*(1+deltaP)
  param1 <- all.parms$model_parms[param]
  
  Y0 = initStates(all.parms$model_parms)
  
  # Run simulation. (Assume original units in mg/L)
  out1 = run_model(mName, times, Y0, all.parms$model_parms)
  out_df1 <- as.data.frame(out1)
  out_df1$time <- out_df1$time/24  # Convert hours to days
  
  #identical(out_df1$time, out_df0$time)
  
  outcome1 <- out_df1$C_bl
  outcome0 <- out_df0$C_bl
  
  #Use AUC
  #outcome1 <- trapz(out_df1$time, out_df1$C_bl)
  #outcome0 <- trapz(out_df0$time, out_df0$C_bl)
  
  #peak concentration
  peakCbl <- out_df0$C_bl[which(out_df0$C_bl == max(out_df0$C_bl))]
  
  #change in outcome:
  ChangeOut <- (outcome1 - outcome0)/outcome0
  #ChangeOut <- (outcome1 - outcome0)/peakCbl
  #change in parameter:
  ChangeParam <- (param1 - param0)/param0
  #sensitivity
  sensitivity <- ChangeOut/ChangeParam
  
  #sensitivity <- ((outcome1 - outcome0)/(all.parms$model_parms[param]*deltaP))*(all.parms$model_parms[param]/outcome0)
  
  #normalize to peak concentration
  
  #sensitivity <- ((outcome1 - outcome0)/(all.parms$model_parms[param]*deltaP))*(all.parms$model_parms[param]/peakCbl)
  
  plot(out_df1$time, sensitivity, type = "l", xlab = "Time [days]", 
       ylab = "Normalized Sensitivity Index", main = param)
  
  plot(out_df1$time, out_df1$C_bl-out_df0$C_bl, type = "l")
  #lines(out_df0$time, out_df0$C_bl, col = "red")
  
  #plot_kim(chem = "PFNA", species = "rat", sex = "male", dose.type = "iv", all.parms, out_df0)
  
  #sensitivity
  #sensitivity[length(sensitivity)]
  
  #peak plasma concentration (returns 0 if using iv dose)
  sensitivity[which(out_df0$C_bl == max(out_df0$C_bl))]
  
}

# set up functions to allow formal fitting of parameters
LS.PBPK.urine <- function(opt.parms, dataVals, dataTimes, mName, all.parms, dose.type, sim.days){
  
  #all.parms["k_ust"] = opt.parms
  
  for (qname in names(opt.parms)) {
    all.parms[qname] = opt.parms[qname]
  }
  
  # Define times for simulation.
  sim_hours = sim.days*24
  #times = seq(from=0.1, to=sim_hours, by=0.1)
  times = dataTimes*24
  #also include times for data points
  #times <- sort(unique(c(times, dataTimes*24)))
  
  all.parms <- initParms(all.parms)
  Y0 = initStates(all.parms)
  
  # Run simulation. (Assume original units in mg/L)
  out = run_model(mName, times, Y0, all.parms)
  out_df <- as.data.frame(out)
  out_df$time <- out_df$time/24  # Convert hours to days
  
  # take only values occuring at dataTimes - this needs checking!!
  #out_dTimes <- out_df[out_df$time %in% dataTimes,]
  modelVals <- (out_df$A_urine/out_df$A_in)*100
  
  if (dose.type == "oral") dataVals <- (dataVals/out_df$A_in)*100
  
  sum((dataVals - modelVals)^2)
  
}

LS.PBPK.all.data <- function(opt.parms,  mName = "PBPK_template", 
                             Udata = NULL, Fdata = NULL, Pdata = NULL,
                             all.parms, log.transform){
  
  if (log.transform == TRUE) opt.parms <- exp(opt.parms)
  
  # Initialize model evaluations terms to NULL
  UmodelVals = NULL
  FmodelVals = NULL
  PmodelVals = NULL
  
  
  for (qname in names(opt.parms)) {
    all.parms[qname] = opt.parms[qname]
  }
  
  # Define times for simulation.  Only use the times from the data (plus zero)
  #times <- sort(unique(c(0, All.Data$Time)))*24
  if (!is.null(Udata)) Utimes <- sort(unique(c(0, Udata[,1])))*24
  if (!is.null(Fdata)) Ftimes <- sort(unique(c(0, Fdata[,1])))*24
  if (!is.null(Pdata)) Ptimes <- sort(unique(c(0, Pdata[,1])))*24
  
  all.parms <- initParms(all.parms)
  Y0 = initStates(all.parms)
  
  # Run simulation. (Assume original units in mg/L)
  # out = run_model(mName, times, Y0, all.parms)
  # out_df <- as.data.frame(out)
  # out_df$time <- out_df$time/24  # Convert hours to days
  
  if (!is.null(Udata)) {
    Uout = run_model(mName, Utimes, Y0, all.parms)
    Uout_df <- as.data.frame(Uout)
    Uout_df$time <- Uout_df$time/24  # Convert hours to days
    Uout_df <- Uout_df[2:length(Uout_df$time),]
    UmodelVals <- (Uout_df$A_urine/Uout_df$A_in)*100
  }
  if (!is.null(Fdata)) {
    Fout = run_model(mName, Ftimes, Y0, all.parms)
    Fout_df <- as.data.frame(Fout)
    Fout_df$time <- Fout_df$time/24  # Convert hours to days
    Fout_df <- Fout_df[2:length(Fout_df$time),]
    FmodelVals <- (Fout_df$A_fecal/Fout_df$A_in)*100
  }
  if (!is.null(Pdata)) {
    Pout = run_model(mName, Ptimes, Y0, all.parms)
    Pout_df <- as.data.frame(Pout)
    Pout_df$time <- Pout_df$time/24  # Convert hours to days
    Pout_df <- Pout_df[2:length(Pout_df$time),]
    PmodelVals <- Pout_df$C_bl  # ug/mL = mg/L
  }
  
  # take only values occuring at dataTimes - this needs checking!! - Unnecessary with current method
  # Uout_dTimes <- out_df[out_df$time %in% UdataVals[,1],]
  # Fout_dTimes <- out_df[out_df$time %in% FdataVals[,1],]
  # Pout_dTimes <- out_df[out_df$time %in% PdataVals[,1],]
  
  # UmodelVals <- (Uout_dTimes$A_urine/Uout_dTimes$A_in)*100
  # FmodelVals <- (Fout_dTimes$A_fecal/Fout_dTimes$A_in)*100
  # #PmodelVals <- Pout_dTimes$C_bl*1000
  
  #if (dose.type == "oral") dataVals <- (dataVals/out_dTimes$A_in)*100
  
  # how to scale data knowing data is sometimes zero? Is scaling necessary? 
  # - might help with pfna female where urine is much higher than feces or plasma
  Usum <- sum((Udata[,2] - UmodelVals)^2)
  Fsum <- sum((Fdata[,2] - FmodelVals)^2)
  Psum <- sum((Pdata[,2] - PmodelVals)^2)
  
  All.sum <- Usum + Fsum + Psum
  All.sum
  
}

LS.PBPK.all.data.mf <- function(opt.parms,  mName = "PBPK_template", 
                                m.Udata = NULL, m.Fdata = NULL, m.Pdata = NULL,
                                f.Udata = NULL, f.Fdata = NULL, f.Pdata = NULL,
                                m.all.parms, f.all.parms, log.transform){
  
  if (log.transform == TRUE) opt.parms <- exp(opt.parms)
  
  # Initialize model evaluations terms to NULL
  f.UmodelVals = NULL
  f.FmodelVals = NULL
  f.PmodelVals = NULL
  
  m.UmodelVals = NULL
  m.FmodelVals = NULL
  m.PmodelVals = NULL
  
  
  # Separate male and female optimization parameters 
  # and substitute in respective all.parms variables
  for (qname in names(opt.parms)) {
    if (str_sub(qname, 1,1) == "f") {
      f.all.parms[str_sub(qname, 3,str_length(qname))] = opt.parms[qname]
    } else if (str_sub(qname, 1,1) == "m") {
      m.all.parms[str_sub(qname, 3,str_length(qname))] = opt.parms[qname]
    } else {
      f.all.parms[qname] = opt.parms[qname]
      m.all.parms[qname] = opt.parms[qname]
    }
  }
  
  # Define times for simulation.  Only use the times from the data (plus zero)
  #times <- sort(unique(c(0, All.Data$Time)))*24
  
  f.times <- sort(unique(c(0, f.Udata[,1], f.Fdata[,1], f.Pdata[,1])))*24
  #if (!is.null(f.Udata)) f.Utimes <- sort(unique(c(0, f.Udata[,1])))*24
  #if (!is.null(f.Fdata)) f.Ftimes <- sort(unique(c(0, f.Fdata[,1])))*24
  #if (!is.null(f.Pdata)) f.Ptimes <- sort(unique(c(0, f.Pdata[,1])))*24
  
  m.times <- sort(unique(c(0, m.Udata[,1], m.Fdata[,1], m.Pdata[,1])))*24
  #if (!is.null(m.Udata)) m.Utimes <- sort(unique(c(0, m.Udata[,1])))*24
  #if (!is.null(m.Fdata)) m.Ftimes <- sort(unique(c(0, m.Fdata[,1])))*24
  #if (!is.null(m.Pdata)) m.Ptimes <- sort(unique(c(0, m.Pdata[,1])))*24
  
  f.all.parms <- initParms(f.all.parms)
  f.Y0 = initStates(f.all.parms)
  
  m.all.parms <- initParms(m.all.parms)
  m.Y0 = initStates(m.all.parms)
  
  
  # Run simulation. (Assume original units in mg/L)
  # out = run_model(mName, times, Y0, all.parms)
  # out_df <- as.data.frame(out)
  # out_df$time <- out_df$time/24  # Convert hours to days
  
  # if (!is.null(f.Udata)) {
  #   f.Uout = run_model(mName, f.Utimes, f.Y0, f.all.parms)
  #   f.Uout_df <- as.data.frame(f.Uout)
  #   f.Uout_df$time <- f.Uout_df$time/24  # Convert hours to days
  #   f.Uout_df <- f.Uout_df[2:length(f.Uout_df$time),]
  #   f.UmodelVals <- (f.Uout_df$A_urine/f.Uout_df$A_in)*100
  # }
  # if (!is.null(f.Fdata)) {
  #   f.Fout = run_model(mName, f.Ftimes, f.Y0, f.all.parms)
  #   f.Fout_df <- as.data.frame(f.Fout)
  #   f.Fout_df$time <- f.Fout_df$time/24  # Convert hours to days
  #   f.Fout_df <- f.Fout_df[2:length(f.Fout_df$time),]
  #   f.FmodelVals <- (f.Fout_df$A_fecal/f.Fout_df$A_in)*100
  # }
  # if (!is.null(f.Pdata)) {
  #   f.Pout = run_model(mName, f.Ptimes, f.Y0, f.all.parms)
  #   f.Pout_df <- as.data.frame(f.Pout)
  #   f.Pout_df$time <- f.Pout_df$time/24  # Convert hours to days
  #   f.Pout_df <- f.Pout_df[2:length(f.Pout_df$time),]
  #   f.PmodelVals <- f.Pout_df$C_bl  # ug/mL = mg/L
  # }
  # if (!is.null(m.Udata)) {
  #   m.Uout = run_model(mName, m.Utimes, m.Y0, m.all.parms)
  #   m.Uout_df <- as.data.frame(m.Uout)
  #   m.Uout_df$time <- m.Uout_df$time/24  # Convert hours to days
  #   m.Uout_df <- m.Uout_df[2:length(m.Uout_df$time),]
  #   m.UmodelVals <- (m.Uout_df$A_urine/m.Uout_df$A_in)*100
  # }
  # if (!is.null(m.Fdata)) {
  #   m.Fout = run_model(mName, m.Ftimes, m.Y0, m.all.parms)
  #   m.Fout_df <- as.data.frame(m.Fout)
  #   m.Fout_df$time <- m.Fout_df$time/24  # Convert hours to days
  #   m.Fout_df <- m.Fout_df[2:length(m.Fout_df$time),]
  #   m.FmodelVals <- (m.Fout_df$A_fecal/m.Fout_df$A_in)*100
  # }
  # if (!is.null(m.Pdata)) {
  #   m.Pout = run_model(mName, m.Ptimes, m.Y0, m.all.parms)
  #   m.Pout_df <- as.data.frame(m.Pout)
  #   m.Pout_df$time <- m.Pout_df$time/24  # Convert hours to days
  #   m.Pout_df <- m.Pout_df[2:length(m.Pout_df$time),]
  #   m.PmodelVals <- m.Pout_df$C_bl  # ug/mL = mg/L
  # }
  
  f.out <- run_model(mName, f.times, f.Y0, f.all.parms)
  f.out.df <- as.data.frame(f.out)
  f.out.df$time <- f.out.df$time/24
  
  f.U.out <- f.out.df[f.out.df$time %in% f.Udata[,1],]
  f.Umodel <- (f.U.out$A_urine/f.U.out$A_in)*100
  f.F.out <- f.out.df[f.out.df$time %in% f.Fdata[,1],]
  f.Fmodel <- (f.F.out$A_fecal/f.F.out$A_in)*100
  f.P.out <- f.out.df[f.out.df$time %in% f.Pdata[,1],]
  f.Pmodel <- f.P.out$C_bl #ug/mL = mg/L
  
  m.out <- run_model(mName, m.times, m.Y0, m.all.parms)
  m.out.df <- as.data.frame(m.out)
  m.out.df$time <- m.out.df$time/24
  
  m.U.out <- m.out.df[m.out.df$time %in% m.Udata[,1],]
  m.Umodel <- (m.U.out$A_urine/m.U.out$A_in)*100
  m.F.out <- m.out.df[m.out.df$time %in% m.Fdata[,1],]
  m.Fmodel <- (m.F.out$A_fecal/m.F.out$A_in)*100
  m.P.out <- m.out.df[m.out.df$time %in% m.Pdata[,1],]
  m.Pmodel <- m.P.out$C_bl #ug/mL = mg/L
  
  f.Usum <- sum((f.Udata[,2] - f.Umodel)^2)
  f.Fsum <- sum((f.Fdata[,2] - f.Fmodel)^2)
  f.Psum <- sum((f.Pdata[,2] - f.Pmodel)^2)
  
  m.Usum <- sum((m.Udata[,2] - m.Umodel)^2)
  m.Fsum <- sum((m.Fdata[,2] - m.Fmodel)^2)
  m.Psum <- sum((m.Pdata[,2] - m.Pmodel)^2)
  
  All.sum <- f.Usum + f.Fsum + f.Psum + m.Usum + m.Fsum + m.Psum
  All.sum
  
}


GLS.PBPK.all.data <- function(opt.parms,  mName = "PBPK_template", 
                              All.data = NULL, all.parms, weights, log.transform){
  
  if (log.transform == TRUE) opt.parms <- exp(opt.parms)
  
  for (qname in names(opt.parms)) {
    all.parms[qname] = opt.parms[qname] 
  }
  
  # Define times for simulation.  Only use the times from the data (plus zero)
  times <- sort(unique(c(0, All.data$Time)))*24
  
  all.parms <- initParms(all.parms)
  Y0 = initStates(all.parms)
  
  out = run_model(mName, times, Y0, all.parms)
  out_df <- as.data.frame(out)
  out_df$time <- out_df$time/24  # Convert hours to days
  out_df <- out_df[2:length(out_df$time),]
  UmodelVals <- (out_df$A_urine/out_df$A_in)*100
  FmodelVals <- (out_df$A_fecal/out_df$A_in)*100
  PmodelVals <- out_df$C_bl  # ug/mL = mg/L
  
  Usum <- sum(weights$Uweights*(All.data$Cum.Urine.Percent - UmodelVals)^2, na.rm = TRUE)
  Fsum <- sum(weights$Fweights*(All.data$Cum.Feces.Percent - FmodelVals)^2, na.rm = TRUE)
  Psum <- sum(weights$Pweights*(All.data$Plasma.Conc - PmodelVals)^2, na.rm = TRUE)
  
  All.sum <- Usum + Fsum + Psum
  All.sum
  
}


PBPK_fit.OLS <- function(mName = "PBPK_template", param.filename = NULL,
                         sheetname = NULL, log.transform = TRUE,
                         opt.parms = c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003)){
  # function to optimize parameters for the PBPK template model
  # Load the dll file
  load_model(mName)
  
  # Adjust default parameters
  all.parms <- load.parameters(mName = "PBPK_template", param.filename = param.filename, sheetname = sheetname)
  
  chem <- all.parms$chem
  species <- all.parms$species
  sex <- all.parms$sex
  dose.type <- all.parms$dose.type
  
  # load urine percent data
  if (dose.type == "oral"){
    dose = all.parms$model_parms["oral_dose_init"]
    AUData <- load_UrineData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    times = AUData[,1]
    U.Data = AUData[,2] #absolute amount, not percent of dose.
  } 
  if (dose.type == "iv"){
    dose = all.parms$model_parms["iv_dose"]
    UData <- load_UrinePercentData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    FData <- load_FecesPercentData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    PData <- load_PlasmaData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    if (all.parms$C.units == "ng/mL") PData[,2] <- PData[,2]/1e3
    PData <- na.omit(PData[,1:2])
  }
  
  #PFNA female params
  # all.parms$model_parms["T_m"] = 0.0100856172
  # all.parms$model_parms["k_loss"] = 0.0001369146
  # all.parms$model_parms["k_ust"] = 0.0294016270
  # all.parms$model_parms["k_bile"] = 1.770978e-04
  # all.parms$model_parms["k_fst"] = 0.036
  
  #PFNA male params
  # all.parms$model_parms["K_t"] = 24.85 #56.6
  # all.parms$model_parms["T_m"] = 7.036256e-03
  # # all.parms$model_parms["k_loss"] = 2.743803e-05
  # all.parms$model_parms["k_bile"] = 2.492374e-05
  # all.parms$model_parms["k_ust"] = 0.00460097
  # all.parms$model_parms["k_u"] = 0.05815641
  
  #PFDA male params
  # all.parms$model_parms["T_m"] = 0.2422122616
  # all.parms$model_parms["k_loss"] = 0.0001310301
  # all.parms$model_parms["k_bile"] = 0.0004101836
  # all.parms$model_parms["k_ust"] = 1.8397075197
  # all.parms$model_parms["k_u"] = 0.007632735
  
  #PFDA female params
  # all.parms$model_parms["T_m"] = 2.040285e-03
  # all.parms$model_parms["k_loss"] = 1.741346e-05
  # all.parms$model_parms["k_bile"] = 2.806921e-04
  # all.parms$model_parms["k_ust"] = 2.383315e-03 #0.002339696 #2.383315e-03
  # all.parms$model_parms["k_u"] = 0.051707393
  
  #opt.parms <- c(T_m = 0.07156389, k_ust = 0.0025)
  #opt.parms <- c(T_m = 0.007957, k_ust = 0.003, k_loss = 0.00003)
  ## Trying to optimize K_t (at the same time as anything other than k_ust) 
  ## results in errors in DLSODA
  
  # opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, F_free = 0.001)
  # opt.parms <- c(T_m = 1, k_ust = 0.001, k_bile = 0.5)#, F_free = 0.001)
  # opt.parms <- c(T_m = runif(1,0,10), k_ust = runif(1,0,10), k_bile = runif(1,0,10))#, F_free = 0.001)
  
  #opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003)
  #opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, K_t = 30.0) 
  # opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003, F_free = 0.00118)
  #Does not work to try to fit everything:
  #opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003, k_u = 1.987, k_fst = 12.56)
  
  #opt.parms <- c(T_m = 0.0079, K_t = 30.0, F_free = 0.00118)
  # opt.parms <- c(k_fst = 2.5e-2, k_u = 0.007632735)
  
  #optimize with Feces data
  # opt.parms <- c(k_fst = 2.5e-2, k_bile = 2.806921e-04)
  #optimize with Urine data
  # opt.parms <- c(k_ust = 2.383315e-03, k_u = 0.00018)
  
  # Optimize multiple parameters just based on urine (percent) data
  #   optResult <- optim(par = opt.parms, fn = LS.PBPK.urine, gr = NULL, method = "Nelder-Mead",
  #                      control=list(parscale = opt.parms),
  # #                     method = "L-BFGS-B", lower = c(0,0,0), upper = c(Inf,Inf,Inf),
  # #                     dataVals = UData[1:6,2], dataTimes = UData[1:6,1], mName = mName,
  #                      dataVals = UData[,2], dataTimes = UData[,1], mName = mName,
  #                      all.parms = all.parms$model_parms, dose.type = dose.type,
  #                      sim.days = all.parms$sim.days)
  
  # Optimize multiple parameters based on urine (percent), feces (percent), and plasma data
  if (log.transform == TRUE) opt.parms = log(opt.parms)
  optResult <- optim(par = opt.parms, fn = LS.PBPK.all.data, gr = NULL, method = "Nelder-Mead",
                     #                     control=list(parscale = opt.parms),
                     #                     method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(Inf,Inf,Inf,Inf),
                     Udata = UData, Fdata = FData, Pdata = PData, mName = mName,
                     # Udata = NULL, Fdata = FData, Pdata = NULL, mName = mName,
                     all.parms = all.parms$model_parms, log.transform = log.transform)
  
  # Optimize single parameter (k_ust) based on urine (percent) data
  #  -requires change to function
  # optResult <- optim(par = 0.83,#, T_m = 34530/1e6, K_t = 24850/1e3),
  #                    fn = LS.PBPK.all.data, gr = NULL, method = "Brent", 
  # lower = 0, upper = 100,
  #                    dataVals = U.Data, dataTimes = times, mName = mName,
  # #                   UdataVals = UData, FdataVals = FData, PdataVals = PData, mName = mName,
  #                    all.parms = all.parms$model_parms, dose.type = dose.type,
  #                    sim.days = all.parms$sim.days)
  
  if (log.transform == TRUE) optResult$exppar <- exp(optResult$par)
  return(optResult)
  
}

PBPK_fit.OLS.mf <- function(mName = "PBPK_template", param.filename = NULL,
                            f.sheetname = NULL, m.sheetname = NULL, log.transform = TRUE,
                            opt.parms = c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003)){
  # function to optimize parameters for the PBPK template model
  # Load the dll file
  load_model(mName)
  
  # Adjust default parameters
  f.all.parms <- load.parameters(mName = mName, param.filename = param.filename, sheetname = f.sheetname)
  m.all.parms <- load.parameters(mName = mName, param.filename = param.filename, sheetname = m.sheetname)
  
  chem <- f.all.parms$chem
  species <- f.all.parms$species
  #sex <- all.parms$sex
  dose.type <- f.all.parms$dose.type
  
  # load data
  f.dose = f.all.parms$model_parms["iv_dose"]
  f.UData <- load_UrinePercentData(chem = chem, species = species, sex = "female", dose.type = dose.type, dose = f.dose)
  f.FData <- load_FecesPercentData(chem = chem, species = species, sex = "female", dose.type = dose.type, dose = f.dose)
  f.PData <- load_PlasmaData(chem = chem, species = species, sex = "female", dose.type = dose.type, dose = f.dose)
  if (f.all.parms$C.units == "ng/mL") f.PData[,2] <- f.PData[,2]/1e3
  f.PData <- na.omit(f.PData[,1:2])
  
  m.dose = m.all.parms$model_parms["iv_dose"]
  m.UData <- load_UrinePercentData(chem = chem, species = species, sex = "male", dose.type = dose.type, dose = m.dose)
  m.FData <- load_FecesPercentData(chem = chem, species = species, sex = "male", dose.type = dose.type, dose = m.dose)
  m.PData <- load_PlasmaData(chem = chem, species = species, sex = "male", dose.type = dose.type, dose = m.dose)
  if (m.all.parms$C.units == "ng/mL") m.PData[,2] <- m.PData[,2]/1e3
  m.PData <- na.omit(m.PData[,1:2])
  
  
  # Optimize multiple parameters based on urine (percent), feces (percent), and plasma data
  if (log.transform == TRUE) opt.parms = log(opt.parms)
  optResult <- optim(par = opt.parms, fn = LS.PBPK.all.data.mf, gr = NULL, method = "Nelder-Mead",
                     #                     control=list(parscale = opt.parms),
                     #                     method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(Inf,Inf,Inf,Inf),
                     f.Udata = f.UData, f.Fdata = f.FData, f.Pdata = f.PData, mName = mName,
                     m.Udata = m.UData, m.Fdata = m.FData, m.Pdata = m.PData,
                     f.all.parms = f.all.parms$model_parms, m.all.parms = m.all.parms$model_parms, log.transform = log.transform)
  
  if (log.transform == TRUE) optResult$exppar <- exp(optResult$par)
  return(optResult)
  
}


PBPK_fit.GLS <- function(mName = "PBPK_template", param.filename = NULL,
                         sheetname = NULL, log.transform = TRUE,
                         opt.parms = c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003)){
  # Function to optimize parameters for the PBPK template model using iterative
  # calling of General Least Squares
  
  # Load the dll file
  load_model(mName)
  
  # Adjust default parameters
  all.parms <- load.parameters(mName = mName, param.filename = param.filename, sheetname = sheetname)
  
  chem <- all.parms$chem
  species <- all.parms$species
  sex <- all.parms$sex
  dose.type <- all.parms$dose.type
  
  # load urine percent data
  if (dose.type == "oral"){
    dose = all.parms$model_parms["oral_dose_init"]
    AUData <- load_UrineData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    times = AUData[,1]
    U.Data = AUData[,2]
  } 
  if (dose.type == "iv"){
    dose = all.parms$model_parms["iv_dose"]
    UData <- load_UrinePercentData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    FData <- load_FecesPercentData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    PData <- load_PlasmaData(chem = chem, species = species, sex = sex, dose.type = dose.type, dose = dose)
    if (all.parms$C.units == "ng/mL") PData[,2] <- PData[,2]/1e3
    PData <- na.omit(PData[,1:2])
    names(PData)[1] <- "Time"
    names(PData)[2] <- "Plasma.Conc"
    All.Data <- merge(x = UData, y = FData, all = TRUE)
    All.Data <- merge(x = All.Data, y = PData, all = TRUE)
  }
  
  #PFNA female params
  # all.parms$model_parms["T_m"] = 0.0100856172
  # all.parms$model_parms["k_loss"] = 0.0001369146
  # all.parms$model_parms["k_ust"] = 0.0294016270
  # all.parms$model_parms["k_bile"] = 1.770978e-04
  # all.parms$model_parms["k_fst"] = 0.036
  #all.parms$model_parms["F_free"] = 0.0015
  
  #PFNA male params
  # all.parms$model_parms["K_t"] = 24.85 #56.6
  # all.parms$model_parms["T_m"] = 7.036256e-03
  # # all.parms$model_parms["k_loss"] = 2.743803e-05
  # all.parms$model_parms["k_bile"] = 2.492374e-05
  # all.parms$model_parms["k_ust"] = 0.00460097
  # all.parms$model_parms["k_u"] = 0.05815641
  
  #PFDA male params
  # all.parms$model_parms["T_m"] = 0.2422122616
  # all.parms$model_parms["k_loss"] = 0.0001310301
  # all.parms$model_parms["k_bile"] = 0.0004101836
  # all.parms$model_parms["k_ust"] = 1.8397075197
  # all.parms$model_parms["k_u"] = 0.007632735
  
  #PFDA female params
  # all.parms$model_parms["T_m"] = 2.040285e-03
  # all.parms$model_parms["k_loss"] = 1.741346e-05
  # all.parms$model_parms["k_bile"] = 2.806921e-04
  # all.parms$model_parms["k_ust"] = 2.383315e-03 #0.002339696 #2.383315e-03
  # all.parms$model_parms["k_u"] = 0.051707393
  
  #opt.parms <- c(T_m = 0.07156389, k_ust = 0.0025)
  #opt.parms <- c(T_m = 0.007957, k_ust = 0.003, k_loss = 0.00003)
  ## Trying to optimize K_t (at the same time as anything other than k_ust) 
  ## results in errors in DLSODA
  
  # opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, F_free = 0.001)
  # opt.parms <- c(T_m = 1, k_ust = 0.001, k_bile = 0.5)#, F_free = 0.001)
  # opt.parms <- c(T_m = runif(1,0,10), k_ust = runif(1,0,10), k_bile = runif(1,0,10))#, F_free = 0.001)
  
  #opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003)
  #opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, K_t = 30.0) 
  # opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003, F_free = 0.00118)
  #Does not work to try to fit everything:
  #opt.parms <- c(T_m = 0.0079, k_ust = 0.003, k_bile = 0.00002, k_loss = 0.00003, k_u = 1.987, k_fst = 12.56)
  
  #opt.parms <- c(T_m = 0.0079, K_t = 30.0, F_free = 0.00118)
  # opt.parms <- c(k_fst = 2.5e-2, k_u = 0.007632735)
  
  #optimize with Feces data
  # opt.parms <- c(k_fst = 2.5e-2, k_bile = 2.806921e-04)
  #optimize with Urine data
  # opt.parms <- c(k_ust = 2.383315e-03, k_u = 0.00018)
  
  # Optimize multiple parameters just based on urine (percent) data
  #   optResult <- optim(par = opt.parms, fn = LS.PBPK.urine, gr = NULL, method = "Nelder-Mead",
  #                      control=list(parscale = opt.parms),
  # #                     method = "L-BFGS-B", lower = c(0,0,0), upper = c(Inf,Inf,Inf),
  # #                     dataVals = UData[1:6,2], dataTimes = UData[1:6,1], mName = mName,
  #                      dataVals = UData[,2], dataTimes = UData[,1], mName = mName,
  #                      all.parms = all.parms$model_parms, dose.type = dose.type,
  #                      sim.days = all.parms$sim.days)
  
  # Step 1: Get initial iterate using ordinary least squares optimization
  # Optimize multiple parameters based on urine (percent), feces (percent), and plasma data
  optResult.old <- optim(par = opt.parms, fn = LS.PBPK.all.data, gr = NULL, method = "Nelder-Mead",
                         # control=list(parscale = opt.parms),
                         # method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(Inf,Inf,Inf,Inf),
                         Udata = UData, Fdata = FData, Pdata = PData, mName = mName,
                         # Udata = NULL, Fdata = FData, Pdata = NULL, mName = mName,
                         all.parms = all.parms$model_parms, log.transform = log.transform)
  
  # do while successive estimates are not sufficiently close
  k <- 1 #iteration counter
  repeat{
    
    # Step 2: Evaluate weights. w = model^2
    for (qname in names(optResult.old$par)) {
      all.parms$model_parms[qname] = optResult.old$par[qname]
    }
    
    times <- sort(unique(c(0, All.Data$Time)))*24
    out_df <- run_sim(mName = mName, parms = all.parms$model_parms, times = times)
    out_df <- out_df[2:length(out_df$time),]
    Uweights <- ((out_df$A_urine/out_df$A_in)*100)^-2
    Fweights <- ((out_df$A_fecal/out_df$A_in)*100)^-2
    Pweights <- (out_df$C_bl)^-2
    weights <- data.frame(All.Data$Time, Uweights, Fweights, Pweights)
    
    # Step 3: Re-estimate q by using weighted least squares optimization
    optResult <- optim(par = opt.parms, fn = GLS.PBPK.all.data, gr = NULL, method = "Nelder-Mead",
                       All.data = All.Data, mName = mName, weights = weights,
                       all.parms = all.parms$model_parms, log.transform = log.transform)
    
    # Check whether successive iterates are close enough
    if ( all( abs(optResult$par - optResult.old$par) < 1e-8 ) ){
      print("Number of iterations:")
      print(k)
      break
    } else if (k > 10000){
      print("Warning: Max number of iterations reached.")
      break
    } else {
      optResult.old <- optResult
      k <- k + 1
    }
    
  }
  
  if (log.transform == TRUE) optResult$exppar <- exp(optResult$par)
  return(optResult)
  
}

