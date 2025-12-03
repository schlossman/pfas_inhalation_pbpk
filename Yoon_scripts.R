# Yoon et al. PBPK model for VOCs (2007)
# Model template simulations 
# Author: Amanda Bernstein, US EPA (ORISE), January 2023

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

# Load functions to run the PBPK model template.
source("run_template_model.R")


# Table 5: Periodic inhalation exposure for 2 month old rats

# VC results
Yoon.table5.VC <- function(){
  # Table 5 VC Results: 
  # 4-hour inhalation exposure, 5 days per week for 2 month old rats
  
  # Set up simulation runs
  conc = c(1, 10000)
  
  Cmax <- rep(0,length(conc)*2)
  AUC <- rep(0,length(conc)*2)
  AM <- rep(0,length(conc)*2)
  AMLVL <- rep(0,length(conc)*2)
  LungEHM <- rep(0,length(conc)*2)
  KidneyEHM <- rep(0,length(conc)*2)
  EHM <- c(rep("no",length(conc)),rep("yes",length(conc)))
  table5 <- data.frame(conc, EHM, Cmax, AUC, AM, AMLVL, LungEHM, KidneyEHM)
  
  idx = 1
  # Compute dose metrics without extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_VC_noEHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_4hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    
    # Compute daily averages over the last week of the simulation
    last <- length(out$C_ven)
    lastweek <- ceiling(length(out$C_ven)/2) #two week simulation
    BW <- 0.187 #kg
    V_li <- BW*0.0492 #kg
    V_li.g <- V_li*1000 #g
    
    table5$Cmax[idx] <- max(out$C_ven)
    table5$AUC[idx] <- (out$AUC_art[last] - out$AUC_art[lastweek])/(7)
    table5$AM[idx] <- (out$A_met_sat[last] - out$A_met_sat[lastweek])/(7)
    table5$AMLVL[idx] <- ((out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7))/V_li.g
    
    idx = idx + 1
  }
  
  # Compute dose metrics with extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_VC_EHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_4hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    
    # Compute daily averages over the last week of the simulation
    last <- length(out$C_ven)
    lastweek <- ceiling(length(out$C_ven)/2) #two week simulation
    BW <- 0.187 #kg
    V_li <- BW*0.0492 #kg
    V_li.g <- V_li*1000 #g
    
    table5$Cmax[idx] <- max(out$C_ven)
    table5$AUC[idx] <- (out$AUC_art[last] - out$AUC_art[lastweek])/(7)
    table5$AM[idx] <- (out$A_met_sat[last] - out$A_met_sat[lastweek])/(7)
    table5$AMLVL[idx] <- ((out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7))/V_li.g
    
    # Extrahepatic metabolism
    livermet <- (out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7)
    lungmet <- (out$A_met_sat_lu[last] - out$A_met_sat_lu[lastweek])/(7)
    kidneymet <- (out$A_met_sat_om[last] - out$A_met_sat_om[lastweek])/(7)
    table5$LungEHM[idx] <- lungmet/livermet*100
    table5$KidneyEHM[idx] <- kidneymet/livermet*100
    
    idx = idx + 1
  }
  
  return(table5)
}

# TCE results
Yoon.table5.TCE <- function(){
  # Table 5 TCE Results: 
  # 8-hour inhalation exposure, 5 days per week for 2 month old rats
  
  # Set up simulation runs
  conc = c(50, 600)
  
  Cmax <- rep(0,length(conc)*2)
  AUC <- rep(0,length(conc)*2)
  AM <- rep(0,length(conc)*2)
  AMLVL <- rep(0,length(conc)*2)
  LungEHM <- rep(0,length(conc)*2)
  KidneyEHM <- rep(0,length(conc)*2)
  EHM <- c(rep("no",length(conc)),rep("yes",length(conc)))
  table5 <- data.frame(conc, EHM, Cmax, AUC, AM, AMLVL, LungEHM, KidneyEHM)
  
  idx = 1
  # Compute dose metrics without extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_TCE_noEHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_8hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    
    # Compute daily averages over the last week of the simulation
    last <- length(out$C_ven)
    lastweek <- ceiling(length(out$C_ven)/2) #four week simulation
    BW <- 0.187 #kg
    V_li <- BW*0.0492 #kg
    V_li.g <- V_li*1000 #g
    
    table5$Cmax[idx] <- max(out$C_ven)
    table5$AUC[idx] <- (out$AUC_art[last] - out$AUC_art[lastweek])/(7)
    table5$AM[idx] <- (out$A_met_sat[last] - out$A_met_sat[lastweek])/(7)
    table5$AMLVL[idx] <- ((out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7))/V_li.g
    
    idx = idx + 1
  }
  
  # Compute dose metrics with extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_TCE_EHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_8hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    
    # Compute daily averages over the last week of the simulation
    last <- length(out$C_ven)
    lastweek <- ceiling(length(out$C_ven)/2) #two week simulation
    BW <- 0.187 #kg
    V_li <- BW*0.0492 #kg
    V_li.g <- V_li*1000 #g
    
    table5$Cmax[idx] <- max(out$C_ven)
    table5$AUC[idx] <- (out$AUC_art[last] - out$AUC_art[lastweek])/(7)
    table5$AM[idx] <- (out$A_met_sat[last] - out$A_met_sat[lastweek])/(7)
    table5$AMLVL[idx] <- ((out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7))/V_li.g
    
    # Extrahepatic metabolism
    livermet <- (out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7)
    lungmet <- (out$A_met_sat_lu[last] - out$A_met_sat_lu[lastweek])/(7)
    kidneymet <- (out$A_met_sat_om[last] - out$A_met_sat_om[lastweek])/(7)
    table5$LungEHM[idx] <- lungmet/livermet*100
    table5$KidneyEHM[idx] <- kidneymet/livermet*100
    
    idx = idx + 1
  }
  
  return(table5)
}

# CCl4 results
Yoon.table5.CCl4 <- function(){
  # Table 5 CCl4 Results: 
  # 6-hour inhalation exposure, 5 days per week for 2 month old rats
  
  # Set up simulation runs
  conc = c(5, 400)
  
  Cmax <- rep(0,length(conc)*2)
  AUC <- rep(0,length(conc)*2)
  AM <- rep(0,length(conc)*2)
  AMLVL <- rep(0,length(conc)*2)
  LungEHM <- rep(0,length(conc)*2)
  KidneyEHM <- rep(0,length(conc)*2)
  EHM <- c(rep("no",length(conc)),rep("yes",length(conc)))
  table5 <- data.frame(conc, EHM, Cmax, AUC, AM, AMLVL, LungEHM, KidneyEHM)
  
  idx = 1
  # Compute dose metrics without extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_CCl4_noEHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_6hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    
    # Compute daily averages over the last week of the simulation
    last <- length(out$C_ven)
    lastweek <- ceiling(length(out$C_ven)/2) #four week simulation
    BW <- 0.187 #kg
    V_li <- BW*0.0492 #kg
    V_li.g <- V_li*1000 #g
    
    table5$Cmax[idx] <- max(out$C_ven)
    table5$AUC[idx] <- (out$AUC_art[last] - out$AUC_art[lastweek])/(7)
    table5$AM[idx] <- (out$A_met_sat[last] - out$A_met_sat[lastweek])/(7)
    table5$AMLVL[idx] <- ((out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7))/V_li.g
    
    idx = idx + 1
  }
  
  # Compute dose metrics with extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_CCl4_EHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_6hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    
    # Compute daily averages over the last week of the simulation
    last <- length(out$C_ven)
    lastweek <- ceiling(length(out$C_ven)/2) #two week simulation
    BW <- 0.187 #kg
    V_li <- BW*0.0492 #kg
    V_li.g <- V_li*1000 #g
    
    table5$Cmax[idx] <- max(out$C_ven)
    table5$AUC[idx] <- (out$AUC_art[last] - out$AUC_art[lastweek])/(7)
    table5$AM[idx] <- (out$A_met_sat[last] - out$A_met_sat[lastweek])/(7)
    table5$AMLVL[idx] <- ((out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7))/V_li.g
    
    # Extrahepatic metabolism
    livermet <- (out$A_met_sat_li[last] - out$A_met_sat_li[lastweek])/(7)
    lungmet <- (out$A_met_sat_lu[last] - out$A_met_sat_lu[lastweek])/(7)
    kidneymet <- (out$A_met_sat_om[last] - out$A_met_sat_om[lastweek])/(7)
    table5$LungEHM[idx] <- lungmet/livermet*100
    table5$KidneyEHM[idx] <- kidneymet/livermet*100
    
    idx = idx + 1
  }
  
  return(table5)
}
