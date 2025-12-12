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
  # As of Dec 8, 2025, the model version with published flow (error in liver
  # blood flow) is not enabled, code commented out, and the corresponding 
  # simulation curves are not created nor shown in the plots.
  
  out <- PBPK_run(model.param.filename = "PFHxS_template_parameters_Model.xlsx", 
                  model.param.sheetname = "FKimRecreateBW", 
                  exposure.param.filename = "PFHxS_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "FKimRecreateBW")

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFHxS.xlsx", 
  #                    sheetname = "FKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc=NULL, chem="PFHxS", sex="female", dose=4.0, img.name=img.name, 
           Pylims=list(plasma=c(0.1,1e5),liver=c(15,30000),kidney=c(1e-4,2e4),urine=c(0,440)))
  
  print.noquote(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print.noquote(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
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
  # As of Dec 8, 2025, the model version with published flow (error in liver
  # blood flow) is not enabled, code commented out, and the corresponding 
  # simulation curves are not created nor shown in the plots.
  
  out <- PBPK_run(model.param.filename = "PFHxS_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MKimRecreateBW", 
                  exposure.param.filename = "PFHxS_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKimRecreateBW",
                  data.times = c(seq(0,1,0.01),seq(1.1,14*24,0.1)))

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFHxS.xlsx", sheetname = "MKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc=NULL, chem="PFHxS", sex="male", dose=10.0, img.name=img.name, 
           Pylims=list(plasma=c(0.1,1.5e5),liver=c(20,50000),kidney=c(1e-4,2e4),urine=c(0,520)))
  
  print.noquote(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print.noquote(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
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
  # As of Dec 8, 2025, the model version with published flow (error in liver
  # blood flow) is not enabled, code commented out, and the corresponding 
  # simulation curves are not created nor shown in the plots.
  
  out <- PBPK_run(model.param.filename = "PFNA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "FKimRecreateBW", 
                  exposure.param.filename = "PFNA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "FKimRecreateBW", 
                  data.times=c(seq(0,1,0.01),seq(1.2,60*24,0.2)))

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFNA.xlsx", sheetname = "FKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc=NULL, chem="PFNA", sex="female", dose=3.0, img.name=img.name, 
           Pylims=list(plasma=c(1e-4,3e1),liver=c(1e-4,20),kidney=c(1e-6,6),urine=c(0,400)))
  
  print.noquote(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print.noquote(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
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
  # As of Dec 8, 2025, the model version with published flow (error in liver
  # blood flow) is not enabled, code commented out, and the corresponding 
  # simulation curves are not created nor shown in the plots.
  
  out <- PBPK_run(model.param.filename = "PFNA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MKimRecreateBW", 
                  exposure.param.filename = "PFNA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKimRecreateBW", 
                  data.times=c(seq(0,1,0.01),seq(1.2,60*24,0.2)))

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFNA.xlsx", sheetname = "MKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc=NULL, chem="PFNA", sex="male", dose=3, img.name=img.name)

  print.noquote(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print.noquote(paste("Maximum mass balance error (Template version with Published Flow):", max(abs(out.inc$A_bal))))
  
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
  # As of Dec 8, 2025, the model version with published flow (error in liver
  # blood flow) is not enabled, code commented out, and the corresponding 
  # simulation curves are not created nor shown in the plots.
  
  out <- PBPK_run(model.param.filename = "PFDA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "FKimRecreateBW", 
                  exposure.param.filename = "PFDA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "FKimRecreateBW", 
                  data.times=c(seq(0,15,0.01),16:(150*24)))

  #out.inc <- PBPK_run(param.filename = "PFAS_template_parameters_PFDA.xlsx", 
  #                   sheetname = "FKimRecreateBW", mName = "PFAS_template_GI")
  plot.Kim(out, out.inc=NULL, chem="PFDA", sex="female", dose=1.0, img.name=img.name)
  
  print.noquote.noquote(paste("Maximum mass balance error (Template version with Corrected Flow):",
              max(abs(out$A_bal))))
  #print.noquote(paste("Maximum mass balance error (Template version with Published Flow):", 
  #     max(abs(out.inc$A_bal))))

  #Accuracy calculation - based on "incorrect" model
  # Pdata <- load_PlasmaData(chem = "PFDA", species = "rat", sex = "female", 
  #                         dose.type = "oral", dose = 1.0)
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

PFOA.Loccisano.KudoLow <- function(img.name=NULL){
  # Figure 8: Loccisano 2012, Male rat, 0.041 mg/kg PFOA (Fig 8 of Loccisano 2012)
  
  out <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MaleRat", 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKudo1BW",
                  data.times = 0:30/10)
  
  plot.Kudo.Loccisano(out, dose="low", img.name=img.name)
  print.noquote(paste("Maximum mass balance error:", max(abs(out$A_bal))))
  
  #Accuracy calculation
  Pdata <- read.csv("Data/Digitized_Data_PFOA/Data_Fig8_Kudo_lowPlasma.csv", 
                    header=TRUE, sep = ",")
  out.inc <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                      model.param.sheetname = "MaleRat", 
                      exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                      exposure.param.sheetname = "MKudo1BW",
                      data.times = Pdata[,3])[-1,] # Without the zero row
  perc <- perc.diff(out.inc$C_bl, Pdata[,4])
  par(mar=c(3, 3, 1, 1), mgp=c(1.5, 0.4, 0))
  plot(out.inc$time, perc, xlab="Time (h)", ylab="% Difference")
  lines(data.times, data.times*0)
  perc.scale <- perc.diff(out.inc$C_bl, Pdata[,4], sc=1.4)
  points(out.inc$time, perc.scale, col="red")
  #within 2.8%
}

PFOA.Loccisano.KudoHigh <- function(img.name = NULL){
  # Figure 9: Loccisano 2012, Male rat, 16.56 mg/kg PFOA
  
  out <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MaleRat", 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKudo2BW")
  
  plot.Kudo.Loccisano(out, dose="high")#, img.name = "Figure9_LoccKudoHigh.tif")
  print.noquote(paste("Maximum mass balance error:", max(abs(out$A_bal))))
  
  #Accuracy calculation
  Pdata <- read.csv("Data/Digitized_Data_PFOA/Data_Fig8_Kudo_highPlasma.csv", 
                    header=TRUE, sep=",")
  nna = !is.na(Pdata[,3])
  out.inc <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                      model.param.sheetname = "MaleRat", 
                      exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                      exposure.param.sheetname = "MKudo2BW",
                      data.times = Pdata[nna,3])[-1,] # Without the zero row
  perc <- perc.diff(out.inc$C_bl, Pdata[nna,4])
  par(mar=c(3, 3, 1, 1), mgp=c(1.5, 0.4, 0))
  plot(out.inc$time, perc, xlab="Time (h)", ylab="% Difference")
  lines(data.times, data.times*0)
  perc.scale <- perc.diff(out.inc$C_bl, Pdata[nna,4], sc=600)
  points(out.inc$time, perc.scale, col="red")
  #within 2.7%
}

PFOA.Loccisano.Kemper <- function(img.name = NULL, match_orig=FALSE){
  # Figure 10: Loccisano 2012, Male rat, 25 mg/kg PFOA
  #
  # The original PBPK Model Template code set the initial condition (IC) for
  # the amount in the GI lumen using the BW in the exposure input spreadsheet,
  # rather than the BW in the exposure-specific BW.table (given below), leading
  # to an initial amount higher than it should  be for this simulation.
  # This error in setting the IC has been corrected in the current code such
  # that the initial BW in the input table is used. However, the resulting 
  # simulation then does not match Figure 10 in Bernstein et al. (2021).
  #
  # *** Set match_orig = TRUE *** as input to thi function applies a dose which
  # repicates the previous dose calculation error, matching that figure.
  #
  # Paul Schlosser, U.S. EPA, Dec 8, 2025
  
  adj.parms = NULL
  if (match_orig)  adj.parms=list(oral_dose_init = 25*0.29/0.206)

  # Construct BW table
  BW.table = list(times = c(0,1,4,7,10,13,17,21)*24, #convert to hours
                  BW = c(0.206,0.264,0.388,0.461,0.510,0.548,0.575,0.590))
  
  out <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                  model.param.sheetname = "MaleRat", 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKemperOral25BW", rtol=1e-14,
                  atol=1e-14, BW.table = BW.table, adj.parms=adj.parms)
  
  plot.Kemper.Loccisano(out, dose.type="oral", dose=25.0, sex="male", img.name=img.name)
  print.noquote(paste("Maximum mass balance error:", max(abs(out$A_bal))))

  #Accuracy calculation
  Pdata <- read.csv("Data/Digitized_Data_PFOA/Data_Fig9_Kemper_OraldosePlasma.csv", 
                    header=TRUE, sep=",")
  out.inc <- PBPK_run(model.param.filename = "PFOA_template_parameters_Model.xlsx", 
                      model.param.sheetname = "MaleRat", 
                      exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                      exposure.param.sheetname = "MKemperOral25BW",
                      BW.table = BW.table, data.times = Pdata[,1])[-1,]
  perc <- perc.diff(out.inc$C_bl, Pdata[,2])
  par(mar=c(3, 3, 1, 1), mgp=c(1.5, 0.4, 0))
  plot(out.inc$time, perc, xlab = "Time (h)", ylab = "% Difference")
  lines(data.times, data.times*0)
  perc.scale <- perc.diff(out.inc$C_bl, Pdata[,2], sc=200)
  points(out.inc$time, perc.scale, col = "red")
  #within 17% near peak, 3% after 70 hours
}

PFOS.Loccisano.3M <- function(img.name = NULL){
  # Figure 11: Loccisano 2012, Male rat, 15 mg/kg PFOS
  
  # Construct BW table
  BW.t =  c( 0, 0.25,  1,  9, 15, 22, 31, 41, 50, 57, 66, 76, 85)*24 # => hours
  BW   = c(233,  233,227,289,330,365,413,456,485,506,527,549,565)/1000 # => kg
  
  # Construct Free fraction table
  times <- c(seq(0,10,0.001),seq(10.001,120*24,0.01)) # 0-120 days => hours
    k_freec = 0.035
    delta = 0.94
    F_free0 = 0.022
  k_freet = k_freec*(approx(BW.t, BW, times)$y^(-0.25))
  Freef.table = list(times = times, 
                     Freef = F_free0*(1 - delta*(1 - exp(-k_freet*times))) )
  
  out <- PBPK_run(model.param.filename = "PFOS_template_parameters_Model.xlsx", 
                  model.param.sheetname = "M3MOralBW", 
                  exposure.param.filename = "PFOS_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "M3MOralBW",
                  BW.table = list(times=BW.t,BW=BW), Freef.table = Freef.table)
  
  plot.3M.Loccisano(out, sex = "male")#, img.name = "Figure11_Locc3MOral2.tif")
  print.noquote(paste("Maximum mass balance error:", max(abs(out$A_bal))))
  
  #Accuracy calculation
  Pdata <- read.csv("Data/Digitized_Data_PFOS/Fig4_3M_OralPlasma.csv", header=TRUE, sep=",")
  out.data <- PBPK_run(model.param.filename = "PFOS_template_parameters_Model.xlsx", 
                       model.param.sheetname = "M3MOralBW", 
                       exposure.param.filename = "PFOS_template_parameters_Exposure.xlsx", 
                       exposure.param.sheetname = "M3MOralBW",
                       BW.table=BW.table, data.times=Pdata[,1]*24, 
                       Freef.table=Freef.table)[]
  out.data <- out.data[-1,] #remove the zero row that was added in PBPK_run()
  perc <- perc.diff(out.data$C_bl, Pdata[,2])
  par(old.par)
  plot(data.times[-(1:2)], perc[-(1:2)], xlab="Time (h)", ylab="% Difference")
  lines(data.times, data.times*0)
  perc <- perc.diff(out.data$C_bl, Pdata[,2], sc=60)
  points(data.times[-(1:2)], perc[-(1:2)], col = "red")
  #within 17%
}
