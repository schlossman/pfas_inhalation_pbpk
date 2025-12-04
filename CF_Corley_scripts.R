#------------------------------------------------------------------------------
# CF_Corley_scripts.R
#
# This file contains functions to recreate results for the Corley chloroform 
# experiments modeled by Sasso et. al. (2013).
#
# Author: Bidya Prasad, December 2021
# Revisions for use with MCSimMod: Pau Schlosser, Nov-Dec 2025
#------------------------------------------------------------------------------

# Set working directory to the directory containing this file.
script.dir=dirname(sys.frame(1)$ofile)
setwd(script.dir)

library(RColorBrewer)
library(readxl)

# Create model and load functions to (compile if needed and) run the PBPK model template.
source("run_template_model.R")

# Corley.chamber.plot contains the commands for the plots in the VOC manuscript submission. 
Corley.chamber.plot<-function(maketiff=FALSE){
  # Creates manuscript Figure 10 (mice) and Figure S-7 (rats)
  # Set input maketiff=TRUE to create tiff file.
  
  vcol=c("#440154FF", "#2A788EFF", "#7AD151FF") 
  vcol5=c("#440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF") 
  vcol7=c("#440154FF", "#472F7DFF", "#39568CFF", "#2A788EFF", "#1F988BFF", 
            "#35B779FF", "#7AD151FF") # viridis(7, begin=0.0, end=0.8)
  
  pub.col <- vcol[2] # also pub.col5.2 <- vcol5[2]
  templ.col <- vcol[3] # also templ.col5.2 <- vcol5[5]
  pub.lty="dashed"; templ.lty="solid"; paper.col=vcol[2]; paper.lty="dashed"
  pub.col5.1=vcol5[1]; pub.col5.2=vcol5[2]; templ.col5.1=vcol5[4];
  templ.col5.2=vcol5[5]; pub.col7.1=vcol7[1]; pub.col7.2=vcol7[2]; pub.col7.3=vcol7[3]; 
  templ.col7.1=vcol7[5]; templ.col7.2=vcol7[6]; templ.col7.3=vcol7[7]
  
  # Simulates 3 rats exposed to chloroform via inhalation in a 9.1L closed 
  # chamber
  conc.ppm <- c(103, 516, 929, 1291, 2581)
  BW <- c(0.226, 0.227, 0.237, 0.239, 0.218)
  # Sasso paper - Corley parameters
  Cor <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Corley_rat_params", 
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname="Corley_rat",
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i]))
    Cor <- cbind(Cor, out$C_chppm)
  }
  Cor <- cbind(Cor, out$time); colnames(Cor) <- c("e1","e2","e3","e4","e5","time")
  
  # Sasso paper - revised parameters
  revi <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Revised_rat_params",
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx",
                    exposure.param.sheetname="Corley_rat",
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i]))
    revi <- cbind(revi, out$C_chppm)
  }
  revi <- cbind(revi, out$time); colnames(revi) <- c("e1","e2","e3","e4","e5","time")
  
  # Load outputs from the acslx Sasso model
  acslx=read_excel("Data/Data_CF/rat_closed_chamber.xlsx")
  # Load inhalation data
  data_rat=read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet="rat")
  
  par(mfrow=c(1,2), mar=c(2.5, 2.5, 1, 1), mgp=c(1.5, 0.4, 0), oma=c(4, 0, 0, 1))
  plot(1,1, type="n", log="y", xlab="Time (hr)", ylab="Chamber Concentration (ppm)", 
       xlim=c(0, 5.5), ylim=c(10, 3000), cex.lab=1.25, yaxt="n", xaxt="n")
  
  axis(side=1, at=seq(0,5.5,0.5), labels=seq(0,5.5,0.5), cex.axis=1)
  axis(side=2, at=c(10,30,100,300,1000,3000), labels=c(10,30,100,300,1000,3000), cex.axis=1)
  points(data_rat$time_hr,data_rat$exp1_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp2_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp3_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp4_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp5_conc_ppm, pch=19, col=pub.col)
  
  # For the E1 scenario  
  lines(Cor[,"time"], Cor[,1], lwd=2, lty=templ.lty, col=templ.col7.1) #1C/T
  lines(acslx$time, acslx$exp1_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #1C/S
  lines(revi[,"time"], revi[,1], lwd=2, lty=templ.lty, col=templ.col) #1R/T
  lines(acslx$time, acslx$exp1_rev, lwd=2, lty=pub.lty, col=pub.col) #1R/S 
  
  # For the E2 scenario
  lines(Cor[,"time"], Cor[,2], lwd=2, lty=templ.lty, col=templ.col7.1) #2C/T
  lines(acslx$time, acslx$exp2_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #2C/S
  lines(revi[,"time"], revi[,2], lwd=2, lty=templ.lty, col=templ.col) #2R/T
  lines(acslx$time, acslx$exp2_rev, lwd=2, lty=pub.lty, col=pub.col) #2R/S 
  
  # For the E3 scenario
  lines(Cor[,"time"], Cor[,3], lwd=2, lty=templ.lty, col=templ.col7.1) #3C/T
  lines(acslx$time, acslx$exp3_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #3C/S
  lines(revi[,"time"], revi[,3], lwd=2, lty=templ.lty, col=templ.col) #3R/T
  lines(acslx$time, acslx$exp3_rev, lwd=2, lty=pub.lty, col=pub.col) #3R/S 
  
  # For the E4 scenario
  lines(Cor[,"time"], Cor[,4], lwd=2, lty=templ.lty, col=templ.col7.1) #4C/T
  lines(acslx$time, acslx$exp4_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #4C/S
  lines(revi[,"time"], revi[,4], lwd=2, lty=templ.lty, col=templ.col) #4R/T
  lines(acslx$time, acslx$exp4_rev, lwd=2, lty=pub.lty, col=pub.col) #4R/S 
  
  # For the E5 scenario
  lines(Cor[,"time"], Cor[,5], lwd=2, lty=templ.lty, col=templ.col7.1) #5C/T
  lines(acslx$time, acslx$exp5_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #5C/S
  lines(revi[,"time"], revi[,5], lwd=2, lty=templ.lty, col=templ.col) #5R/T
  lines(acslx$time, acslx$exp5_rev, lwd=2, lty=pub.lty, col=pub.col) #5R/S 
  
  
  # Simulates 15 mice exposed to chloroform via inhalation in a 9.1L closed 
  # chamber
  conc.ppm <- c(1000, 2500, 5000); BW <- c(0.028, 0.029, 0.029)
  # Sasso paper - Corley parameters
  Cor <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Corley_mouse_params",
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx",
                    exposure.param.sheetname="Corley_mouse",
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i]))
    Cor <- cbind(Cor, out$C_chppm)
  }
  Cor <- cbind(Cor, out$time); colnames(Cor) <- c("e1","e2","e3", "time")
  
  # Sasso paper - revised parameters
  revi <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Revised_mouse_params", 
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname="Corley_mouse",
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i]))
    revi <- cbind(revi, out$C_chppm)
  }
  revi <- cbind(revi, out$time); colnames(revi) <- c("e1","e2","e3","time")
  
  # Load outputs from the acslx Sasso model
  acslx=read_excel("Data/Data_CF/mouse_closed_chamber.xlsx")
  # Load inhalation data
  data_mouse=read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet="mouse")
  
  plot(1,1, type="n", log="y", xlab="Time (hr)",cex.lab=1.25, yaxt="n", xaxt="n",  
       ylab="Chamber Concentration (ppm)", xlim=c(0, 5.5), ylim=c(10, 5000))
  
  axis(side=1, at=seq(0,5.5,0.5), labels=seq(0,5.5,0.5), cex.axis=1)
  axis(side=2, at=c(10,30,100,300,1000,3000), labels=c(10,30,100,300,1000,3000), cex.axis=1)
  points(data_mouse$time_hr,data_mouse$exp1_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp2_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp3_conc_ppm, pch=19, col=pub.col)
  
  # For the E1 scenario  
  lines(Cor[,"time"], Cor[,1], lwd=2, lty=templ.lty, col=templ.col7.1) #1C/T
  lines(acslx$time, acslx$exp1_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #1C/S
  lines(revi[,"time"], revi[,1], lwd=2, lty=templ.lty, col=templ.col) #1R/T
  lines(acslx$time, acslx$exp1_rev, lwd=2, lty="dotdash", col=pub.col) #1R/S 
  
  # For the E2 scenario
  lines(Cor[,"time"], Cor[,2], lwd=2, lty=templ.lty, col=templ.col7.1) #2C/T
  lines(acslx$time, acslx$exp2_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #2C/S
  lines(revi[,"time"], revi[,2], lwd=2, lty=templ.lty, col=templ.col) #2R/T
  lines(acslx$time, acslx$exp2_rev, lwd=2, lty="dotdash", col=pub.col) #2R/S 
  
  # For the E3 scenario
  lines(Cor[,"time"], Cor[,3], lwd=2, lty=templ.lty, col=templ.col7.1) #3C/T
  lines(acslx$time, acslx$exp3_Cor, lwd=2, lty=pub.lty, col=pub.col5.2) #3C/S
  lines(revi[,"time"], revi[,3], lwd=2, lty=templ.lty, col=templ.col) #3R/T
  lines(acslx$time, acslx$exp3_rev, lwd=2, lty="dotdash", col=pub.col) #3R/S 
  
  par(mfrow=c(1,1), oma=c(0, 0, 1, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
  legend("bottomleft", c("Published Data", "Corley assumptions (Published Model)", 
                         "Corley assumptions (Template Version)"),
         xpd=TRUE, inset=c(0.05,0), bty="n", col=c(pub.col,pub.col5.2,templ.col7.1),
         lty=c(NA,pub.lty,templ.lty), pch=c(19,NA,NA), lwd=c(2,2,2), cex=1.12)
  legend("bottomleft", c("Revised assumptions (Published Model)", 
                         "Revised assumptions (Template Version)"),
         xpd=TRUE, inset=c(0.55, 0), bty="n", col=c(pub.col, templ.col),
         lty=c("dotdash", templ.lty), pch=c(NA,NA), lwd=c(2,2), cex=1.12)
  
  # Plot mouse simulation: Figure 10
  if (maketiff){
  tiff("CF_Corley_mouse.tiff", res=300, height=8, width=8, units="in")
  par(mar=c(5, 5, 2, 0), oma=c(4, 0.5, 0.5, 2.5))
  plot(1,1, type="n", log="y", xlab="Time (hr)", cex.lab=1.25, yaxt="n", xaxt="n", 
       ylab="Chamber Concentration (ppm)", xlim=c(0,5.5), ylim=c(10, 5000))
  
  axis(side=1, at=seq(0,5.5,0.5), labels=seq(0,5.5,0.5), cex.axis=1)
  axis(side=2, at=c(10,100,1000), labels=c(10,100,1000), cex.axis=1)
  points(data_mouse$time_hr,data_mouse$exp1_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp2_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp3_conc_ppm, pch=19, col=pub.col)
  
  # For the E1 scenario  
  lines(Cor[,"time"], Cor[,1], lwd=4, lty=templ.lty, col=templ.col7.1) #1C/T
  lines(acslx$time, acslx$exp1_Cor, lwd=4, lty=pub.lty, col=pub.col5.2) #1C/S
  lines(revi[,"time"], revi[,1], lwd=4, lty=templ.lty, col=templ.col) #1R/T
  lines(acslx$time, acslx$exp1_rev, lwd=4, lty="dotdash", col=pub.col) #1R/S 
  
  # For the E2 scenario
  lines(Cor[,"time"], Cor[,2], lwd=4, lty=templ.lty, col=templ.col7.1) #2C/T
  lines(acslx$time, acslx$exp2_Cor, lwd=4, lty=pub.lty, col=pub.col5.2) #2C/S
  lines(revi[,"time"], revi[,2], lwd=4, lty=templ.lty, col=templ.col) #2R/T
  lines(acslx$time, acslx$exp2_rev, lwd=4, lty="dotdash", col=pub.col) #2R/S 
  
  # For the E3 scenario
  lines(Cor[,"time"], Cor[,3], lwd=4, lty=templ.lty, col=templ.col7.1) #3C/T
  lines(acslx$time, acslx$exp3_Cor, lwd=4, lty=pub.lty, col=pub.col5.2) #3C/S
  lines(revi[,"time"], revi[,3], lwd=4, lty=templ.lty, col=templ.col) #3R/T
  lines(acslx$time, acslx$exp3_rev, lwd=4, lty="dotdash", col=pub.col) #3R/S 
  
  par(mfrow=c(1,1), oma=c(0, 0, 1, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
  legend("bottomleft", c("Published Data", "Corley assumptions (Published Model)", 
                         "Corley assumptions (Template Version)"),
         xpd=TRUE, inset=c(0.08, 0.02), bty="n", pch=c(19,NA,NA), lwd=c(2,2,3),
         lty=c(NA, pub.lty, templ.lty), col=c(pub.col, pub.col5.2, templ.col7.1))
  legend("bottomleft", c("Revised assumptions (Published Model)", 
                         "Revised assumptions (Template Version)"),
         xpd=TRUE, inset=c(0.52, 0.02), bty="n", pch=c(NA,NA), lwd=c(2,3),
         col=c(pub.col, templ.col), lty=c("dotdash", templ.lty))
  dev.off()
  } # End if (maketiff)
}

# plot.Corley.pub contains the commands for the plots in the VOC manuscript submission.
# This function should be run after the CF.Corley.chamber.total.plot.pub() function.
# The plot is saved as a vector file (".svg") in the user's working directory.

plot.Corley.pub <- function(){
  svg(filename="CF_Corley_erat_mouse.svg", width=12, height=5, pointsize=12)
  CF.Corley.chamber.total.plot.pub()
  dev.off()
}

Corley.rat.diffs <- function(){
  # Calculate error: percent difference between the template and Sasso 
  # acslx models for: (1) Corley parameters and (2) revised parameters
  # Discard first data point at time=0, C_ven=0
  
  # Import acslx results 
  acslx <- as.data.frame(read_excel("Data/Data_CF/rat_closed_chamber.xlsx"))
  
  conc.ppm <- c(103, 516, 929, 1291, 2581); BW <- c(0.226, 0.227, 0.237, 0.239, 0.218)
   
  # Sasso paper - Corley parameters
  Cor <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Corley_rat_params", 
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname="Corley_rat",
                    data.times=acslx$time,
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i]))
    Cor <- cbind(Cor, out$C_chppm)
  }
  Cor <- as.data.frame(cbind(Cor, out$time)); colnames(Cor) <- c("e1","e2","e3","e4","e5","time")
  
  # Sasso paper - revised parameters
  revi <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Revised_rat_params",
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx",
                    exposure.param.sheetname="Corley_rat",
                    data.times=acslx$time,
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i]))
    revi <- cbind(revi, out$C_chppm)
  }
  revi <- as.data.frame(cbind(revi, out$time)); colnames(revi) <- c("e1","e2","e3","e4","e5","time")

  Cv=match(acslx$time,Cor$time); rv=match(acslx$time,revi$time)
  print("Maximum percentage difference between Template and acslx rat models",quote=FALSE)
  print("Experiment 1",quote=FALSE)
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,1], acslx$exp1_Cor)),quote=FALSE)
  print(paste0("b. Revised parameters: ", max.diff(revi[rv,1], acslx$exp1_rev)),quote=FALSE)
  print("Experiment 2",quote=FALSE)
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,2], acslx$exp2_Cor)),quote=FALSE)
  print(paste("b. Revised parameters:", max.diff(revi[rv,2], acslx$exp2_rev)),quote=FALSE)
  print("Experiment 3")
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,3], acslx$exp3_Cor)),quote=FALSE)
  print(paste("b. Revised parameters:", max.diff(revi[rv,3], acslx$exp3_rev)),quote=FALSE)
  print("Experiment 4",quote=FALSE)
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,4], acslx$exp4_Cor)),quote=FALSE)
  print(paste("b. Revised parameters:", max.diff(revi[rv,4], acslx$exp4_rev)),quote=FALSE)
  print("Experiment 5",quote=FALSE)
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,5], acslx$exp5_Cor)),quote=FALSE)
  print(paste("b. Revised parameters:", max.diff(revi[rv,5], acslx$exp5_rev)),quote=FALSE)
}

Corley.mouse.diffs <- function(){
  # Calculate error: percent difference between the template and Sasso 
  # acslx models for: (1) Corley parameters and (2) revised parameters
  
  # Import acslX resuts 
  acslx <- as.data.frame(read_excel("Data/Data_CF/mouse_closed_chamber.xlsx"))
  
  conc.ppm <- c(1000, 2500, 5000); BW <- c(0.028, 0.029, 0.029)
  # Sasso paper - Corley parameters
  Cor <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Corley_mouse_params", 
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname="Corley_mouse",
                    data.times=acslx[["time"]],
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i], NCH=15))
    Cor <- cbind(Cor, out$C_chppm)
  }
  Cor <- as.data.frame(cbind(Cor,out$time)); colnames(Cor) <- c("e1","e2","e3","time")
  
  # Sasso paper - revised parameters
  revi <- NULL
  for (i in 1:length(conc.ppm)){
    out <- PBPK_run(model.param.filename="CF_template_parameters_Model.xlsx",
                    model.param.sheetname="Revised_mouse_params", 
                    exposure.param.filename="CF_template_parameters_Exposure.xlsx", 
                    exposure.param.sheetname="Corley_mouse",
                    data.times=acslx[["time"]],
                    adj.parms=c(Conc_init=conc.ppm[i], BW=BW[i], NCH=15))
    revi <- cbind(revi, out$C_chppm)
  }
  revi <- as.data.frame(cbind(revi, out$time)); colnames(revi) <- c("e1","e2","e3","time")

  Cv=match(acslx$time,Cor$time); rv=match(acslx$time,revi$time)
  print("Maximum percentage difference between Template and acslx mouse models",quote=FALSE)
  print("Experiment 1",quote=FALSE)
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,1], acslx$exp1_Cor)),quote=FALSE)
  print(paste0("b. Revised parameters: ", max.diff(revi[rv,1], acslx$exp1_rev)),quote=FALSE)
  print("Experiment 2",quote=FALSE)
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,2], acslx$exp2_Cor)),quote=FALSE)
  print(paste("b. Revised parameters:", max.diff(revi[rv,2], acslx$exp2_rev)),quote=FALSE)
  print("Experiment 3")
  print(paste("a. Corley parameters:", max.diff(Cor[Cv,3], acslx$exp3_Cor)),quote=FALSE)
  print(paste("b. Revised parameters:", max.diff(revi[rv,3], acslx$exp3_rev)),quote=FALSE)
}