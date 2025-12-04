# Yoon et al. PBPK model for VOCs (2007)
# Model template simulations 
# Author: Amanda Bernstein, US EPA (ORISE), January 2023
# Revisions for use with MCSimMod: Paul Schlosser, December 2025

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile); setwd(script.dir)
source("run_template_model.R") # Load functions to run the PBPK model template.

# Table 5: Periodic inhalation exposure for 2 month old rats

make_t5 <- function(conc){ # Creates table 5 (t5) structure filled with NaNs
  Cmax <- AUC <- AM <- AMLVL <- LungEHM <- KidneyEHM <- rep(NaN,length(conc))
  EHM <- c(rep("no",length(conc)),rep("yes",length(conc)))
  return(data.frame(conc, EHM, Cmax, AUC, AM, AMLVL, LungEHM, KidneyEHM))
}

daily_avg <- function(out, t5, idx){
  # Compute daily averages over the last week of simulations in Yoon.table5 
  # functions below
  nt <- length(out$time)
  nl <- max(which(out$time.days <= (out$time.days[nt]-7))) # Start of last week
  t5$Cmax[idx] <- max(out$C_ven)
  t5$AUC[idx] <- (out$AUC_art[nt] - out$AUC_art[nl])/7
  t5$AM[idx] <- (out$A_met_sat[nt] - out$A_met_sat[nl])/7
  t5$AMLVL[idx] <- (out$A_met_sat_li[nt] - out$A_met_sat_li[nl])/(7*out$V_li[nt]*1000)
      # AMLVL is normalized to liver mass in grams
  # Extra-hepatic metabolism
  livermet <- out$A_met_sat_li[nt] - out$A_met_sat_li[nl]
  t5$LungEHM[idx] <- 100*(out$A_met_sat_lu[nt] - out$A_met_sat_lu[nl])/livermet
  t5$KidneyEHM[idx] <- 100*(out$A_met_sat_om[nt] - out$A_met_sat_om[nl])/livermet
  return(t5)
}

# VC results
Yoon.table5.VC <- function(){
  # Table 5 VC Results: 
  # 4-hour inhalation exposure, 5 days per week for 2 month old rats
  
  # Load values saved from previous version of this function (used the old MCSim method):
  s5 <- as.data.frame(readRDS("Data/Yoon.table5.VC.rds"))
  
  # Set up simulation runs
  t5 <- make_t5(conc <- c(1,10000)); idx = 1
  
  # Compute dose metrics without extra-hepatic metabolism
  for (ii in conc){ # for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_VC_noEHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_4hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    t5 <- daily_avg(out, t5, idx)
    idx = idx + 1
  }
  
  # Compute dose metrics with extra-hepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_VC_EHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_4hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    t5 <- daily_avg(out, t5, idx)
    idx = idx + 1
  }
  # Append table of % differences between current t5 and saved version (s5),
  # rounded to three significant figures:
  print("Replication of Yoon et al. (2007) Table 5 values (first 4 rows) for VC",quote=FALSE)
  print("and % differences between previously saved calculations (last 4 rows).",quote=FALSE)
  s5[,3:8] <- round(100*(s5[,3:8] - t5[,3:8])/(s5[,3:8]+1e-18), 3)
  t5[,3:8] <- signif(t5[,3:8],3)
  return(rbind(t5,s5))
}

# TCE results
Yoon.table5.TCE <- function(){
  # Table 5 TCE Results: 
  # 8-hour inhalation exposure, 5 days per week for 2 month old rats
  
  # Load values saved from previous version of this function (used the old MCSim method):
  s5 <- as.data.frame(readRDS("Data/Yoon.table5.TCE.rds"))
  
  # Set up simulation runs
  t5 <- make_t5(conc <- c(50,600)); idx = 1
  
  # Compute dose metrics without extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_TCE_noEHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_8hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    t5 <- daily_avg(out, t5, idx)
    idx = idx + 1
  }
  
  # Compute dose metrics with extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_TCE_EHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_8hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    t5 <- daily_avg(out, t5, idx)
    idx = idx + 1
  }
  # Append table of % differences between current t5 and saved version (s5),
  # rounded to three significant figuress:
  print("Replication of Yoon et al. (2007) Table 5 values (first 4 rows) for TCE",quote=FALSE)
  print("and % differences between previously saved calculations (last 4 rows).",quote=FALSE)
  s5[,3:8] <- round(100*(s5[,3:8] - t5[,3:8])/(s5[,3:8]+1e-18), 3)
  t5[,3:8] <- signif(t5[,3:8],3)
  return(rbind(t5,s5))
}

# CCl4 results
Yoon.table5.CCl4 <- function(){
  # Table 5 CCl4 Results: 
  # 6-hour inhalation exposure, 5 days per week for 2 month old rats
  
  # Load values saved from previous version of this function (used the old MCSim method):
  s5 <- as.data.frame(readRDS("Data/Yoon.table5.CCl4.rds"))
  
  # Set up simulation runs
  t5 <- make_t5(conc <- c(5,400)); idx = 1
  
  # Compute dose metrics without extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_CCl4_noEHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_6hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    t5 <- daily_avg(out, t5, idx)
    idx = idx + 1
  }
  
  # Compute dose metrics with extrahepatic metabolism
  for (ii in conc){ #for each concentration
    out <- PBPK_run(model.param.filename = "Yoon_template_parameters_Model.xlsx",
                    model.param.sheetname = "2mo_CCl4_EHM", 
                    exposure.param.filename = "Yoon_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname = "inhal_6hr_5dweek", 
                    adj.parms = c(Conc_init = ii))
    t5 <- daily_avg(out, t5, idx)
    idx = idx + 1
  }
  # Append table of % differences between current t5 and saved version (s5),
  # rounded to three significant figuress:
  print("Replication of Yoon et al. (2007) Table 5 values (first 4 rows) for CCl4",quote=FALSE)
  print("and % differences between previously saved calculations (last 4 rows).",quote=FALSE)
  s5[,3:8] <- round(100*(s5[,3:8] - t5[,3:8])/(s5[,3:8]+1e-18), 3)
  t5[,3:8] <- signif(t5[,3:8],3)
  return(rbind(t5,s5))
  
}
