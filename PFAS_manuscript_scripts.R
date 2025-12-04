# Cases run for PFAS PBPK template manuscript:
#  Bernstein, et al. (2021)
# Requires source of "run_template_model.R" and "plot_PFAS_template_man.R"
#
# Note, this file has been updated to use the updated version of the PBPK 
# template. For the models of PFHxS, PFNA, and PFDA from Kim et al. (2018) and 
# (2019), this file does not include simulations using the incorrect published 
# flow rates leaving the liver. Therefore, those simulation lines do not appear
# in the figures and the accuracy calculations which rely on those values have 
# been suppressed. To see the simulations with the incorrect published flow 
# rate leaving the liver for these models, see the PFAS PBPK template published
# on the US EPA Environmental Dataset Gateway (https://doi.org/10.23719/1520081).

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

old.par <- par(no.readonly = TRUE)
source("run_template_model.R")
source("plot_PFAS_template_man.R")

PFHxS.Kim.FemaleRat <- function(img.name = NULL){
  # Figure 3: Kim 2018, Female rat, 4 mg/kg PFHxS
  
  out <- PBPK_run(model.param.filename = "PFHxS_template_parameters_Model.xlsx", 
                  model.param.sheetname = "FKimRecreateBW", 
                  exposure.param.filename = "PFHxS_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "FKimRecreateBW")

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFHxS.xlsx", sheetname = "FKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc = NULL, chem = "PFHxS", sex = "female", dose = 4.0, img.name = img.name)
  
  print(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
  #Accuracy calculation - based on "incorrect" model
  # Pdata <- load_PlasmaData(chem = "PFHxS", species = "rat", sex = "female", dose.type = "oral", dose = 4.0)
  # data.times <- Pdata[,3]*24
  # out.inc.data <- PBPK_run(model.param.filename = "PFHxS_template_parameters_Model.xlsx", 
  #                          model.param.sheetname = "FKimRecreateBW", 
  #                          exposure.param.filename = "PFHxS_template_parameters_Exposure.xlsx", 
  #                          exposure.param.sheetname = "FKimRecreateBW",
  #                          mName = "PFAS_template_GI", 
  #                          data.times = data.times)
  # out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  # perc <- perc.diff(model = out.inc.data$C_bl*1e3, data = Pdata[,4])
  # par(old.par)
  # plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference", ylim = c(-0.1,5.5))
  # lines(data.times, data.times*0)
  # perc.scale <- perc.diff.scale(model = out.inc.data$C_bl*1e3, data = Pdata[,4], fig.scale = 10^6-1)
  # points(data.times, perc.scale, col = "red")
  # #within 0.1%
}

PFHxS.Kim.MaleRat <- function(img.name = NULL){
  # Figure 4: Kim 2018, Male rat, 10 mg/kg PFHxS
  
  out <- PBPK_run(model.param.filename = "PFHxS_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MKimRecreateBW", 
                  exposure.param.filename = "PFHxS_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKimRecreateBW")

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFHxS.xlsx", sheetname = "MKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc = NULL, chem = "PFHxS", sex = "male", dose = 10.0, img.name = img.name)
  
  print(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
  #Accuracy calculation - based on "incorrect" model
  # Pdata <- load_PlasmaData(chem = "PFHxS", species = "rat", sex = "male", dose.type = "oral", dose = 10.0)
  # data.times <- Pdata[,3]*24
  # out.inc.data <- PBPK_run(model.param.filename = "PFHxS_template_parameters_Model.xlsx", 
  #                          model.param.sheetname = "MKimRecreateBW", 
  #                          exposure.param.filename = "PFHxS_template_parameters_Exposure.xlsx", 
  #                          exposure.param.sheetname = "MKimRecreateBW",
  #                          mName = "PFAS_template_GI", 
  #                          data.times = data.times)
  # out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  # perc <- perc.diff(model = out.inc.data$C_bl*1e3, data = Pdata[,4])
  # par(old.par)
  # plot(data.times[3:length(data.times)], perc[3:length(data.times)], type = "p", xlab = "Time (h)", ylab = "% Difference")
  # lines(data.times, data.times*0)
  # perc.scale <- perc.diff.scale(model = out.inc.data$C_bl*1e3, data = Pdata[,4], fig.scale = 10^6-10^2)
  # points(data.times, perc.scale, col = "red")
  # #within 0.3%
}

PFNA.Kim.FemaleRat <- function(img.name = NULL){
  # Figure 5: Kim 2019, Female rat, 3 mg/kg PFNA
  
  out <- PBPK_run(model.param.filename = "PFNA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "FKimRecreateBW", 
                  exposure.param.filename = "PFNA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "FKimRecreateBW")

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFNA.xlsx", sheetname = "FKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc = NULL, chem = "PFNA", sex = "female", dose = 3.0, img.name = img.name)
  
  print(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
  #Accuracy calculation - based on "incorrect" model
  # Pdata <- load_PlasmaData(chem = "PFNA", species = "rat", sex = "female", dose.type = "oral", dose = 3.0)
  # data.times <- Pdata[,3]*24
  # out.inc.data <- PBPK_run(model.param.filename = "PFNA_template_parameters_Model.xlsx", 
  #                          model.param.sheetname = "FKimRecreateBW", 
  #                          exposure.param.filename = "PFNA_template_parameters_Exposure.xlsx", 
  #                          exposure.param.sheetname = "FKimRecreateBW",
  #                          mName = "PFAS_template_GI", 
  #                          data.times = data.times)
  # out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  # perc <- perc.diff(model = out.inc.data$C_bl, data = Pdata[,4])
  # par(old.par)
  # plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference")
  # lines(data.times, data.times*0)
  # perc.scale <- perc.diff.scale(model = out.inc.data$C_bl, data = Pdata[,4], fig.scale = 100-0.0001)
  # points(data.times, perc.scale, col = "red")
  # #within 0.1%
}

PFNA.Kim.MaleRat <- function(img.name = NULL){
  # Figure 6: Kim 2019, Male rat, 3 mg/kg PFNA
  out <- PBPK_run(model.param.filename = "PFNA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MKimRecreateBW", 
                  exposure.param.filename = "PFNA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKimRecreateBW")

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFNA.xlsx", sheetname = "MKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc = NULL, chem = "PFNA", sex = "male", dose = 3.0, img.name = img.name)

  print(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
  #Accuracy calculation - based on "incorrect" model
  # Pdata <- load_PlasmaData(chem = "PFNA", species = "rat", sex = "male", dose.type = "oral", dose = 3.0)
  # data.times <- Pdata[,3]*24
  # out.inc.data <- PBPK_run(model.param.filename = "PFNA_template_parameters_Model.xlsx", 
  #                          model.param.sheetname = "MKimRecreateBW", 
  #                          exposure.param.filename = "PFNA_template_parameters_Exposure.xlsx", 
  #                          exposure.param.sheetname = "MKimRecreateBW",
  #                          mName = "PFAS_template_GI", 
  #                          data.times = data.times)
  # out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  # perc <- perc.diff(model = out.inc.data$C_bl, data = Pdata[,4]) #warning due to NA values in data
  # par(old.par)
  # plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference")
  # lines(data.times, data.times*0)
  # perc.scale <- perc.diff.scale(model = out.inc.data$C_bl, data = Pdata[,4], fig.scale = 100-0.0001)
  # points(data.times, perc.scale, col = "red")
  # #within 0.3%
}

PFDA.Kim.FemaleRat <- function(img.name = NULL){
  # Figure 7: Kim 2019, Female rat, 1 mg/kg PFDA
  
  out <- PBPK_run(model.param.filename = "PFDA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "FKimRecreateBW", 
                  exposure.param.filename = "PFDA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "FKimRecreateBW")

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFDA.xlsx", sheetname = "FKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc = NULL, chem = "PFDA", sex = "female", dose = 1.0, img.name = img.name)
  
  print(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
  #Accuracy calculation - based on "incorrect" model
  # Pdata <- load_PlasmaData(chem = "PFDA", species = "rat", sex = "female", dose.type = "oral", dose = 1.0)
  # data.times <- Pdata[,3]*24
  # out.inc.data <- PBPK_run(model.param.filename = "PFDA_template_parameters_Model.xlsx", 
  #                          model.param.sheetname = NULL, 
  #                          exposure.param.filename = "PFDA_template_parameters_Exposure.xlsx", 
  #                          exposure.param.sheetname = NULL,
  #                          #sheetname = "FKimRecreateBW", mName = "PFAS_template_GI", 
  #                          data.times = data.times)
  # out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  # perc <- perc.diff(model = out.inc.data$C_bl, data = Pdata[,4])
  # par(old.par)
  # plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference")
  # lines(data.times, data.times*0)
  # perc.scale <- perc.diff.scale(model = out.inc.data$C_bl, data = Pdata[,4], fig.scale = 10-0.01)
  # points(data.times, perc.scale, col = "red")
  # #within 0.6%
}

PFOA.Loccisano.KudoLow <- function(img.name = NULL){
  # Figure 8: Loccisano 2012, Male rat, 0.041 mg/kg PFOA (Fig 8 of Loccisano 2012)
  
  out <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MaleRat", 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKudo1BW")
  
  plot.Kudo.Loccisano(out, dose = "low", img.name = img.name)
  print(paste("Maximum mass balance error:", max(abs(out$A_bal))))
  
  #Accuracy calculation
  data.loc <- "Data/Digitized_Data_PFOA/"
  Pdata <- read.csv(file = paste0(data.loc, "Data_Fig8_Kudo_lowPlasma.csv"), header = TRUE, sep = ",")
  data.times <- Pdata[,3]
  out.inc.data <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                           model.param.sheetname = "MaleRat", 
                           exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                           exposure.param.sheetname = "MKudo1BW",
                           data.times = data.times)
  out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  perc <- perc.diff(model = out.inc.data$C_bl, data = Pdata[,4])
  par(old.par)
  plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference")
  lines(data.times, data.times*0)
  perc.scale <- perc.diff.scale(model = out.inc.data$C_bl, data = Pdata[,4], fig.scale = 1.4)
  points(data.times, perc.scale, col = "red")
  #within 2.8%
}

PFOA.Loccisano.KudoHigh <- function(img.name = NULL){
  # Figure 9: Loccisano 2012, Male rat, 16.56 mg/kg PFOA
  
  out <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MaleRat", 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKudo2BW")
  
  plot.Kudo.Loccisano(out, dose = "high")#, img.name = "Figure9_LoccKudoHigh.tif")
  print(paste("Maximum mass balance error:", max(abs(out$A_bal))))
  
  #Accuracy calculation
  data.loc <- "Data/Digitized_Data_PFOA/"
  Pdata <- read.csv(file = paste0(data.loc, "Data_Fig8_Kudo_highPlasma.csv"), header = TRUE, sep = ",")
  data.times <- Pdata[,3]
  out.inc.data <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                           model.param.sheetname = "MaleRat", 
                           exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                           exposure.param.sheetname = "MKudo2BW",
                           data.times = data.times)
  out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  perc <- perc.diff(model = out.inc.data$C_bl, data = Pdata[,4])
  par(old.par)
  plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference")
  lines(data.times, data.times*0)
  perc.scale <- perc.diff.scale(model = out.inc.data$C_bl, data = Pdata[,4], fig.scale = 600)
  points(data.times, perc.scale, col = "red")
  #within 2.7%
}

PFOA.Loccisano.Kemper <- function(img.name = NULL){
  # Figure 10: Loccisano 2012, Male rat, 25 mg/kg PFOA
  
  # Construct BW table
  times = c(0,1,4,7,10,13,17,21)*24 #convert to hours
  BW = c(0.206,0.264,0.388,0.461,0.510,0.548,0.575,0.590)
  BW.table = list(times = times, BW = BW)
  
  out <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MaleRat", 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKemperOral25BW",
                  BW.table = BW.table)
  
  plot.Kemper.Loccisano(out, dose.type = "oral", dose = 25.0, sex = "male", img.name = img.name)
  print(paste("Maximum mass balance error:", max(abs(out$A_bal))))

  #Accuracy calculation
  data.loc <- "Data/Digitized_Data_PFOA/"
  Pdata <- read.csv(file = paste0(data.loc, "Data_Fig9_Kemper_OraldosePlasma.csv"), header = TRUE, sep = ",")
  data.times <- Pdata[,1]
  out.inc.data <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                           model.param.sheetname = "MaleRat", 
                           exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                           exposure.param.sheetname = "MKemperOral25BW",
                           BW.table = BW.table, data.times = data.times)
  out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  perc <- perc.diff(model = out.inc.data$C_bl, data = Pdata[,2])
  par(old.par)
  plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference")
  lines(data.times, data.times*0)
  perc.scale <- perc.diff.scale(model = out.inc.data$C_bl, data = Pdata[,2], fig.scale = 200)
  points(data.times, perc.scale, col = "red")
  #within 17% near peak, 3% after 70 hours
}

PFOS.Loccisano.3M <- function(img.name = NULL){
  # Figure 11: Loccisano 2012, Male rat, 15 mg/kg PFOS
  
  # Construct BW table
  BW_times = c(0, 0.25, 1, 9,  15,  22, 31,  41, 50, 57, 66, 76, 85)*24 #convert to hours
  BW = c(0.233,0.233,0.227,0.289,0.330,0.365,0.413,0.456,0.485,0.506,0.527,0.549,0.565)
  BW.table = list(times = BW_times, BW = BW)
  
  # Construct Free fraction table
  endtime = 120*24 # convert 120 days to hours
  Freetimes1 = seq(from=0, to=10, by=0.00001)
  Freetimes2 = seq(from=10, to=endtime, by=0.001)
  Freetimes = sort(unique(c(Freetimes1, Freetimes2)))
  BW_fnc <- approxfun(BW_times, BW)
  BW_freetimes <- BW_fnc(Freetimes)
  k_freec = 0.035
  delta = 0.94
  F_free = 0.022
  k_free = k_freec*(BW_freetimes^-0.25)
  Free <- F_free*(1 - delta*(1 - exp(-k_free*Freetimes)))
  Freef.table = list(times = Freetimes, Freef = Free)
  
  out <- PBPK_run(model.param.filename = "PFOS_template_parameters_Model.xlsx", 
                  model.param.sheetname = "M3MOralBW", 
                  exposure.param.filename = "PFOS_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "M3MOralBW",
                  BW.table = BW.table, Freef.table = Freef.table)
  
  plot.3M.Loccisano(out, sex = "male")#, img.name = "Figure11_Locc3MOral2.tif")
  print(paste("Maximum mass balance error:", max(abs(out$A_bal))))
  
  #Accuracy calculation
  data.loc <- "Data/Digitized_Data_PFOS/"
  Pdata <- read.csv(file = paste0(data.loc, "Fig4_3M_OralPlasma.csv"), header = TRUE, sep = ",")
  data.times <- Pdata[,1]*24
  out.data <- PBPK_run(model.param.filename = "PFOS_template_parameters_Model.xlsx", 
                       model.param.sheetname = "M3MOralBW", 
                       exposure.param.filename = "PFOS_template_parameters_Exposure.xlsx", 
                       exposure.param.sheetname = "M3MOralBW",
                       BW.table = BW.table, data.times = data.times, Freef.table = Freef.table)
  out.data <- out.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  perc <- perc.diff(model = out.data$C_bl, data = Pdata[,2])
  par(old.par)
  plot(data.times[3:length(data.times)], perc[3:length(data.times)], type = "p", xlab = "Time (h)", ylab = "% Difference")
  lines(data.times, data.times*0)
  perc.scale <- perc.diff.scale(model = out.data$C_bl, data = Pdata[,2], fig.scale = 60)
  points(data.times[3:length(data.times)], perc.scale[3:length(data.times)], col = "red")
  #within 1%
}