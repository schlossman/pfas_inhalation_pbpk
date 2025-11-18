#------------------------------------------------------------------------------
# CF_Corley_scripts.R
#
# This file contains functions to recreate results for the Corley chloroform 
# experiments modeled by Sasso et. al. (2013).
#
# Author: Bidya Prasad, December 2021
# Revisions for use with MCSimMod: Pau Schlosser, November 2025
#------------------------------------------------------------------------------

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

library(RColorBrewer)
library(readxl)

# Create model and load functions to (compile if needed and) run the PBPK model template.
source("run_template_model.R")

# CF.Corley.chamber.total.plot.pub contains the commands for the plots in the VOC
# manuscript submission. 
CF.Corley.chamber.total.plot.pub <-function(){
  # Creates manuscript Figure 10 (mice) and Figure S-7 (rats)
  
  vcol = c("#440154FF", "#2A788EFF", "#7AD151FF") 
  vcol5 = c("#440154FF", "#414487FF", "#2A788EFF", "#22A884FF", "#7AD151FF") 
  vcol7 = c("#440154FF", "#472F7DFF", "#39568CFF", "#2A788EFF", "#1F988BFF", 
            "#35B779FF", "#7AD151FF") # viridis(7, begin = 0.0, end = 0.8)
  
  pub.col <- vcol[2] # also pub.col5.2 <- vcol5[2]
  templ.col <- vcol[3] # also templ.col5.2 <- vcol5[5]
  pub.lty <- "dashed"
  templ.lty <- "solid"
  paper.col <- vcol[2]
  paper.lty <- "dashed"
  
  pub.col5.1 <- vcol5[1]
  pub.col5.2 <- vcol5[2]
  templ.col5.1 <- vcol5[4]
  templ.col5.2 <- vcol5[5]
  
  pub.col7.1 <- vcol7[1]
  pub.col7.2 <- vcol7[2]
  pub.col7.3 <- vcol7[3]
  templ.col7.1 <- vcol7[5]
  templ.col7.2 <- vcol7[6]
  templ.col7.3 <- vcol7[7]
  
  
  # Simulates 3 rats exposed to chloroform via inhalation in a 9.1L closed 
  # chamber
  conc.ppm <- c(103, 516, 929, 1291, 2581)
  BW <- c(0.226, 0.227, 0.237, 0.239, 0.218)
  # Sasso paper - Corley parameters
  out.all.Cor.param <- NULL
  for (i in 1:nd){
    out.Cor.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Corley_rat_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_rat",
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i]))
    out.all.Cor.param <- cbind(out.all.Cor.param, out.Cor.param$C_chppm)
  }
  colnames(out.all.Cor.param) <- c("exp_1","exp_2","exp_3","exp_4","exp_5")
  
  # Sasso paper - revised parameters
  out.all.rev.param <- NULL
  for (i in 1:length(conc.ppm)){
    out.rev.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Revised_rat_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_rat",
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i]))
    out.all.rev.param <- cbind(out.all.rev.param, out.rev.param$C_chppm)
  }
  colnames(out.all.rev.param) <- c("exp_1","exp_2","exp_3","exp_4","exp_5")
  
  # Load outputs from the acslx Sasso model
  acslx_rat_out = read_excel("Data/Data_CF/rat_closed_chamber.xlsx")
  # Load inhalation data
  data_rat = read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet = "rat")
  
  par(mfrow=c(1,2), mar=c(2.5, 2.5, 2.0, 1.1), mgp = c(1.5, 0.4, 0), oma = c(3.2, 0, 0, 1))
  plot(1,1, type="n", log="y", 
       xlab="Time (hr)",  
       ylab="Chamber Concentration (ppm)", xlim=c(0, 5.5), ylim=c(10, 5000), 
       cex.lab = 1.25, yaxt="n", xaxt="n")
  
  axis(side=1, at=seq(from=0, to=5.5, by=0.5), 
       labels=c(0, 0.5, 1, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 5.5), 
       cex.axis = 1)
  axis(side=2, at=c(1, 10, 100, 1000), labels=c(1, 10, 100, 1000), 
       cex.axis = 1)
  points(data_rat$time_hr,data_rat$exp1_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp2_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp3_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp4_conc_ppm, pch=19, col=pub.col)
  points(data_rat$time_hr,data_rat$exp5_conc_ppm, pch=19, col=pub.col)
  
  # For the Exp_1 scenario  
  lines(out.Cor.param$time, out.all.Cor.param[,1], lwd = 2, lty=templ.lty, col=templ.col7.1) #1C/T
  lines(acslx_rat_out$time, acslx_rat_out$exp1_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #1C/S
  lines(out.rev.param$time, out.all.rev.param[,1], lwd = 2, lty=templ.lty, col=templ.col) #1R/T
  lines(acslx_rat_out$time, acslx_rat_out$exp1_rev, lwd = 2, lty=pub.lty, col=pub.col) #1R/S 
  
  # For the Exp_2 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,2], lwd = 2, lty=templ.lty, col=templ.col7.1) #2C/T
  lines(acslx_rat_out$time, acslx_rat_out$exp2_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #2C/S
  lines(out.rev.param$time, out.all.rev.param[,2], lwd = 2, lty=templ.lty, col=templ.col) #2R/T
  lines(acslx_rat_out$time, acslx_rat_out$exp2_rev, lwd = 2, lty=pub.lty, col=pub.col) #2R/S 
  
  # For the Exp_3 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,3], lwd = 2, lty=templ.lty, col=templ.col7.1) #3C/T
  lines(acslx_rat_out$time, acslx_rat_out$exp3_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #3C/S
  lines(out.rev.param$time, out.all.rev.param[,3], lwd = 2, lty=templ.lty, col=templ.col) #3R/T
  lines(acslx_rat_out$time, acslx_rat_out$exp3_rev, lwd = 2, lty=pub.lty, col=pub.col) #3R/S 
  
  # For the Exp_4 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,4], lwd = 2, lty=templ.lty, col=templ.col7.1) #4C/T
  lines(acslx_rat_out$time, acslx_rat_out$exp4_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #4C/S
  lines(out.rev.param$time, out.all.rev.param[,4], lwd = 2, lty=templ.lty, col=templ.col) #4R/T
  lines(acslx_rat_out$time, acslx_rat_out$exp4_rev, lwd = 2, lty=pub.lty, col=pub.col) #4R/S 
  
  # For the Exp_5 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,5], lwd = 2, lty=templ.lty, col=templ.col7.1) #5C/T
  lines(acslx_rat_out$time, acslx_rat_out$exp5_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #5C/S
  lines(out.rev.param$time, out.all.rev.param[,5], lwd = 2, lty=templ.lty, col=templ.col) #5R/T
  lines(acslx_rat_out$time, acslx_rat_out$exp5_rev, lwd = 2, lty=pub.lty, col=pub.col) #5R/S 
  
  
  # Simulates 15 mice exposed to chloroform via inhalation in a 9.1L closed 
  # chamber
  conc.ppm <- c(1000, 2500, 5000)
  BW <- c(0.028, 0.029, 0.029)
  # Sasso paper - Corley parameters
  out.all.Cor.param <- NULL
  for (i in 1:length(conc.ppm)){
    out.Cor.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Corley_mouse_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_mouse",
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i]))
    out.all.Cor.param <- cbind(out.all.Cor.param, out.Cor.param$C_chppm)
  }
  colnames(out.all.Cor.param) <- c("exp_1","exp_2","exp_3")
  
  # Sasso paper - revised parameters
  out.all.rev.param <- NULL
  for (i in 1:length(conc.ppm)){
    out.rev.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Revised_mouse_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_mouse",
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i]))
    out.all.rev.param <- cbind(out.all.rev.param, out.rev.param$C_chppm)
  }
  colnames(out.all.rev.param) <- c("exp_1","exp_2","exp_3")
  
  # Load outputs from the acslx Sasso model
  acslx_mouse_out = read_excel("Data/Data_CF/mouse_closed_chamber.xlsx")
  # Load inhalation data
  data_mouse = read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet = "mouse")
  
  plot(1,1, type="n", log="y", 
       xlab="Time (hr)",  
       ylab="Chamber Concentration (ppm)", xlim=c(0, 5.5), ylim=c(10, 5000), 
       cex.lab = 1.25, yaxt="n", xaxt="n")
  
  axis(side=1, at=seq(from=0, to=5.5, by=0.5), 
       labels=c(0, 0.5, 1, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 5.5), 
       cex.axis = 1)
  axis(side=2, at=c(10, 100, 1000), labels=c(10, 100, 1000), 
       cex.axis = 1)
  points(data_mouse$time_hr,data_mouse$exp1_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp2_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp3_conc_ppm, pch=19, col=pub.col)
  
  # For the Exp_1 scenario  
  lines(out.Cor.param$time, out.all.Cor.param[,1], lwd = 2, lty=templ.lty, col=templ.col7.1) #1C/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp1_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #1C/S
  lines(out.rev.param$time, out.all.rev.param[,1], lwd = 2, lty=templ.lty, col=templ.col) #1R/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp1_rev, lwd = 2, lty="dotdash", col=pub.col) #1R/S 
  
  # For the Exp_2 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,2], lwd = 2, lty=templ.lty, col=templ.col7.1) #2C/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp2_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #2C/S
  lines(out.rev.param$time, out.all.rev.param[,2], lwd = 2, lty=templ.lty, col=templ.col) #2R/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp2_rev, lwd = 2, lty="dotdash", col=pub.col) #2R/S 
  
  # For the Exp_3 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,3], lwd = 2, lty=templ.lty, col=templ.col7.1) #3C/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp3_Cor, lwd = 2, lty=pub.lty, col=pub.col5.2) #3C/S
  lines(out.rev.param$time, out.all.rev.param[,3], lwd = 2, lty=templ.lty, col=templ.col) #3R/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp3_rev, lwd = 2, lty="dotdash", col=pub.col) #3R/S 
  
  par(mfrow = c(1,1), oma = c(0, 0, 1, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottomleft", c("Published Data", "Corley assumptions (Published Model)", "Corley assumptions (Template Version)"),
         xpd = TRUE, inset = c(0.1, -0.025), bty = "n",
         col = c(pub.col, pub.col5.2, templ.col7.1),
         lty = c(NA, pub.lty, templ.lty),
         pch = c(19,NA,NA), lwd = c(2,2,2), cex=1.12)
  legend("bottomleft", c("Revised assumptions (Published Model)", "Revised assumptions (Template Version)"),
         xpd = TRUE, inset = c(0.57, -0.01), bty = "n",
         col = c(pub.col, templ.col),
         lty = c("dotdash", templ.lty),
         pch = c(NA,NA), lwd = c(2,2), cex=1.12)
  
  # Plot mouse simulation: Figure 10
  tiff("CF_Corley_mouse.tiff", res=300, height=8, width=8, units="in")
  par(mar = c(5, 5, 2, 0), oma = c(4, 0.5, 0.5, 2.5))
  plot(1,1, type="n", log="y", 
       xlab="Time (hr)",  
       ylab="Chamber Concentration (ppm)", xlim=c(0, 5.5), ylim=c(10, 5000), 
       cex.lab = 1.25, yaxt="n", xaxt="n")
  
  axis(side=1, at=seq(from=0, to=5.5, by=0.5), 
       labels=c(0, 0.5, 1, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5, 5.5), 
       cex.axis = 1)
  axis(side=2, at=c(10, 100, 1000), labels=c(10, 100, 1000), 
       cex.axis = 1)
  points(data_mouse$time_hr,data_mouse$exp1_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp2_conc_ppm, pch=19, col=pub.col)
  points(data_mouse$time_hr,data_mouse$exp3_conc_ppm, pch=19, col=pub.col)
  
  # For the Exp_1 scenario  
  lines(out.Cor.param$time, out.all.Cor.param[,1], lwd = 4, lty=templ.lty, col=templ.col7.1) #1C/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp1_Cor, lwd = 4, lty=pub.lty, col=pub.col5.2) #1C/S
  lines(out.rev.param$time, out.all.rev.param[,1], lwd = 4, lty=templ.lty, col=templ.col) #1R/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp1_rev, lwd = 4, lty="dotdash", col=pub.col) #1R/S 
  
  # For the Exp_2 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,2], lwd = 4, lty=templ.lty, col=templ.col7.1) #2C/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp2_Cor, lwd = 4, lty=pub.lty, col=pub.col5.2) #2C/S
  lines(out.rev.param$time, out.all.rev.param[,2], lwd = 4, lty=templ.lty, col=templ.col) #2R/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp2_rev, lwd = 4, lty="dotdash", col=pub.col) #2R/S 
  
  # For the Exp_3 scenario
  lines(out.Cor.param$time, out.all.Cor.param[,3], lwd = 4, lty=templ.lty, col=templ.col7.1) #3C/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp3_Cor, lwd = 4, lty=pub.lty, col=pub.col5.2) #3C/S
  lines(out.rev.param$time, out.all.rev.param[,3], lwd = 4, lty=templ.lty, col=templ.col) #3R/T
  lines(acslx_mouse_out$time, acslx_mouse_out$exp3_rev, lwd = 4, lty="dotdash", col=pub.col) #3R/S 
  
  par(mfrow = c(1,1), oma = c(0, 0, 1, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottomleft", c("Published Data", "Corley assumptions (Published Model)", "Corley assumptions (Template Version)"),
         xpd = TRUE, inset = c(0.08, 0.02), bty = "n",
         col = c(pub.col, pub.col5.2, templ.col7.1),
         lty = c(NA, pub.lty, templ.lty),
         pch = c(19,NA,NA), lwd = c(2,2,3))
  legend("bottomleft", c("Revised assumptions (Published Model)", "Revised assumptions (Template Version)"),
         xpd = TRUE, inset = c(0.52, 0.02), bty = "n",
         col = c(pub.col, templ.col),
         lty = c("dotdash", templ.lty),
         pch = c(NA,NA), lwd = c(2,3))
  
  dev.off()
}

# plot.Corley.pub contains the commands for the plots in the VOC manuscript submission.
# This function should be run after the CF.Corley.chamber.total.plot.pub() function.
# The plot is saved as a vector file (".svg") in the user's working directory.

plot.Corley.pub <- function(){
  svg(filename="CF_Corley_exp_rat_mouse.svg", 
    width=12, 
    height=5, 
    pointsize=12)
CF.Corley.chamber.total.plot.pub()
dev.off()
}

CF.Corley.percent.rat.diffs <- function(){
  # Calculate error: percent difference between the template and Sasso 
  # acslx models for: (1) Corley parameters and (2) revised parameters
  # Discard first data point at time = 0, C_ven = 0
  
  # Import data 
  data_rat = read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet = "rat")
  acslx_rat_out = read_excel("Data/Data_CF/rat_closed_chamber.xlsx")
  
  
  conc.ppm <- c(103, 516, 929, 1291, 2581)
  BW <- c(0.226, 0.227, 0.237, 0.239, 0.218)
   
  # Sasso paper - Corley parameters
  out.all.Cor.param <- NULL
  for (i in 1:length(conc.ppm)){
    out.Cor.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Corley_rat_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_rat",
                              data.times = acslx_rat_out[["time"]],
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i]))
    out.all.Cor.param <- cbind(out.all.Cor.param, out.Cor.param$C_chppm)
  }
  colnames(out.all.Cor.param) <- c("exp_1","exp_2","exp_3","exp_4","exp_5")
  
  # Sasso paper - revised parameters
  out.all.rev.param <- NULL
  for (i in 1:length(conc.ppm)){
    out.rev.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Revised_rat_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_rat",
                              data.times = acslx_rat_out[["time"]],
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i]))
    out.all.rev.param <- cbind(out.all.rev.param, out.rev.param$C_chppm)
  }
  colnames(out.all.rev.param) <- c("exp_1","exp_2","exp_3","exp_4","exp_5")

  data_rat = read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet = "rat")
  acslx_rat_out = read_excel("Data/Data_CF/rat_closed_chamber.xlsx")
  
  ntimes_acslx = length(acslx_rat_out[["time"]])
  ntimes_template = length(out.rev.param$time)

  df_C_exp_1_temp <- data.frame(time = out.rev.param$time,
                                exp_1_Cor_temp = out.all.Cor.param[,1],
                                exp_1_rev_temp = out.all.rev.param[,1])
  df_C_exp_1_acslx <- data.frame(time = acslx_rat_out["time"],
                                 exp_1_Cor_acslx = acslx_rat_out["exp1_Cor"],
                                 acslx_rat_out["exp1_rev"])
  df_acslx_temp_merge_exp_1 <- merge(df_C_exp_1_temp, df_C_exp_1_acslx, by.x='time')
  
  df_C_exp_2_temp = data.frame(time = out.rev.param$time,
                               exp_2_Cor_temp = out.all.Cor.param[,2],
                               exp_2_rev_temp = out.all.rev.param[,2])
  df_C_exp_2_acslx <- data.frame(time = acslx_rat_out["time"],
                                 exp_2_Cor_acslx = acslx_rat_out["exp2_Cor"],
                                 acslx_rat_out["exp2_rev"])
  df_acslx_temp_merge_exp_2 <- merge(df_C_exp_2_temp, df_C_exp_2_acslx, by.x='time')
  
  df_C_exp_3_temp = data.frame(time = out.rev.param$time,
                               exp_3_Cor_temp = out.all.Cor.param[,3],
                               exp_3_rev_temp = out.all.rev.param[,3])
  df_C_exp_3_acslx <- data.frame(time = acslx_rat_out["time"],
                                 exp_3_Cor_acslx = acslx_rat_out["exp3_Cor"],
                                 acslx_rat_out["exp3_rev"])
  df_acslx_temp_merge_exp_3 <- merge(df_C_exp_3_temp, df_C_exp_3_acslx, by.x='time')
  
  df_C_exp_4_temp = data.frame(time = out.rev.param$time,
                               exp_4_Cor_temp = out.all.Cor.param[,4],
                               exp_4_rev_temp = out.all.rev.param[,4])
  df_C_exp_4_acslx <- data.frame(time = acslx_rat_out["time"],
                                 exp_4_Cor_acslx = acslx_rat_out["exp4_Cor"],
                                 acslx_rat_out["exp4_rev"])
  df_acslx_temp_merge_exp_4 <- merge(df_C_exp_4_temp, df_C_exp_4_acslx, by.x='time')
  
  df_C_exp_5_temp = data.frame(time = out.rev.param$time,
                               exp_5_Cor_temp = out.all.Cor.param[,5],
                               exp_5_rev_temp = out.all.rev.param[,5])
  df_C_exp_5_acslx <- data.frame(time = acslx_rat_out["time"],
                                 exp_5_Cor_acslx = acslx_rat_out["exp5_Cor"],
                                 acslx_rat_out["exp5_rev"])
  df_acslx_temp_merge_exp_5 <- merge(df_C_exp_5_temp, df_C_exp_5_acslx, by.x='time')
  

  err_exp_1_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_1[,"exp_1_Cor_temp"], 
                                               df_acslx_temp_merge_exp_1[,"exp1_Cor"])))
  err_exp_1_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_1[,"exp_1_rev_temp"], 
                                               df_acslx_temp_merge_exp_1[,"exp1_rev"])))
  
  print("Maximum percentage difference between Template and acslx models")
  print("1. Experiment 1")
  print(paste0("a.    Corley parameters: ", err_exp_1_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_1_rev_temp_acslx))
  
  err_exp_2_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_2[,"exp_2_Cor_temp"], 
                                               df_acslx_temp_merge_exp_2[,"exp2_Cor"])))
  err_exp_2_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_2[,"exp_2_rev_temp"], 
                                               df_acslx_temp_merge_exp_2[,"exp2_rev"])))
  print("2. Experiment 2")
  print(paste0("a.    Corley parameters: ", err_exp_2_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_2_rev_temp_acslx))
  
  err_exp_3_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_3[,"exp_3_Cor_temp"], 
                                               df_acslx_temp_merge_exp_3[,"exp3_Cor"])))
  err_exp_3_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_3[,"exp_3_rev_temp"], 
                                               df_acslx_temp_merge_exp_3[,"exp3_rev"])))
  print("3. Experiment 3")
  print(paste0("a.    Corley parameters: ", err_exp_3_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_3_rev_temp_acslx))
  
  err_exp_4_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_4[,"exp_4_Cor_temp"], 
                                               df_acslx_temp_merge_exp_4[,"exp4_Cor"])))
  err_exp_4_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_4[,"exp_4_rev_temp"], 
                                               df_acslx_temp_merge_exp_4[,"exp4_rev"])))
  print("4. Experiment 4")
  print(paste0("a.    Corley parameters: ", err_exp_4_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_4_rev_temp_acslx))
  
  err_exp_5_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_5[,"exp_5_Cor_temp"], 
                                               df_acslx_temp_merge_exp_5[,"exp5_Cor"])))
  err_exp_5_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_5[,"exp_5_rev_temp"], 
                                               df_acslx_temp_merge_exp_5[,"exp5_rev"])))
  print("5. Experiment 5")
  print(paste0("a.    Corley parameters: ", err_exp_5_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_5_rev_temp_acslx))
}

CF.Corley.percent.mouse.diffs <- function(){
  # Calculate error: percent difference between the template and Sasso 
  # acslx models for: (1) Corley parameters and (2) revised parameters
  
  # Import data 
  data_mouse = read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet = "mouse")
  acslx_mouse_out = read_excel("Data/Data_CF/mouse_closed_chamber.xlsx")
  
  conc.ppm <- c(1000, 2500, 5000)
  BW <- c(0.028, 0.029, 0.029)
  # Sasso paper - Corley parameters
  out.all.Cor.param <- NULL
  for (i in 1:length(conc.ppm)){
    out.Cor.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Corley_mouse_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_mouse",
                              data.times = acslx_mouse_out[["time"]],
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i], NCH = 15))
    out.all.Cor.param <- cbind(out.all.Cor.param, out.Cor.param$C_chppm)
  }
  colnames(out.all.Cor.param) <- c("exp_1","exp_2","exp_3")
  
  # Sasso paper - revised parameters
  out.all.rev.param <- NULL
  for (i in 1:length(conc.ppm)){
    out.rev.param <- PBPK_run(model.param.filename = "CF_template_parameters_Model.xlsx",
                              model.param.sheetname = "Revised_mouse_params", 
                              exposure.param.filename = "CF_template_parameters_Exposure.xlsx", 
                              exposure.param.sheetname = "Corley_mouse",
                              data.times = acslx_mouse_out[["time"]],
                              adj.parms = c(Conc_init = conc.ppm[i], BW = BW[i], NCH = 15))
    out.all.rev.param <- cbind(out.all.rev.param, out.rev.param$C_chppm)
  }
  colnames(out.all.rev.param) <- c("exp_1","exp_2","exp_3")

  data_mouse = read_excel("Data/Data_CF/corley_1990_chamber_exp.xlsx", sheet = "mouse")
  acslx_mouse_out = read_excel("Data/Data_CF/mouse_closed_chamber.xlsx")
  ntimes_acslx = length(acslx_mouse_out[["time"]])
  ntimes_template = length(out.rev.param$time)
  
  df_C_exp_1_temp <- data.frame(time = out.rev.param$time,
                                exp_1_Cor_temp = out.all.Cor.param[,1],
                                exp_1_rev_temp = out.all.rev.param[,1])
  df_C_exp_1_acslx <- data.frame(time = acslx_mouse_out["time"],
                                 exp_1_Cor_acslx = acslx_mouse_out["exp1_Cor"],
                                 acslx_mouse_out["exp1_rev"])
  df_acslx_temp_merge_exp_1 <- merge(df_C_exp_1_temp, df_C_exp_1_acslx, by.x='time')
  
  df_C_exp_2_temp = data.frame(time = out.rev.param$time,
                               exp_2_Cor_temp = out.all.Cor.param[,2],
                               exp_2_rev_temp = out.all.rev.param[,2])
  df_C_exp_2_acslx <- data.frame(time = acslx_mouse_out["time"],
                                 exp_2_Cor_acslx = acslx_mouse_out["exp2_Cor"],
                                 acslx_mouse_out["exp2_rev"])
  df_acslx_temp_merge_exp_2 <- merge(df_C_exp_2_temp, df_C_exp_2_acslx, by.x='time')
  
  df_C_exp_3_temp = data.frame(time = out.rev.param$time,
                               exp_3_Cor_temp = out.all.Cor.param[,3],
                               exp_3_rev_temp = out.all.rev.param[,3])
  df_C_exp_3_acslx <- data.frame(time = acslx_mouse_out["time"],
                                 exp_3_Cor_acslx = acslx_mouse_out["exp3_Cor"],
                                 acslx_mouse_out["exp3_rev"])
  df_acslx_temp_merge_exp_3 <- merge(df_C_exp_3_temp, df_C_exp_3_acslx, by.x='time')
  
  err_exp_1_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_1[,"exp_1_Cor_temp"], 
                                               df_acslx_temp_merge_exp_1[,"exp1_Cor"])))
  err_exp_1_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_1[,"exp_1_rev_temp"], 
                                               df_acslx_temp_merge_exp_1[,"exp1_rev"])))
  print("Maximum percentage difference between Template and acslx models")
  print("1. Experiment 1")
  print(paste0("a.    Corley parameters: ", err_exp_1_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_1_rev_temp_acslx))
  
  err_exp_2_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_2[,"exp_2_Cor_temp"], 
                                               df_acslx_temp_merge_exp_2[,"exp2_Cor"])))
  err_exp_2_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_2[,"exp_2_rev_temp"], 
                                               df_acslx_temp_merge_exp_2[,"exp2_rev"])))
  print("2. Experiment 2")
  print(paste0("a.    Corley parameters: ", err_exp_2_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_2_rev_temp_acslx))
  
  err_exp_3_Cor_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_3[,"exp_3_Cor_temp"], 
                                               df_acslx_temp_merge_exp_3[,"exp3_Cor"])))
  err_exp_3_rev_temp_acslx = max(abs(perc.diff(df_acslx_temp_merge_exp_3[,"exp_3_rev_temp"], 
                                               df_acslx_temp_merge_exp_3[,"exp3_rev"])))
  print("3. Experiment 3")
  print(paste0("a.    Corley parameters: ", err_exp_3_Cor_temp_acslx))
  print(paste0("b.    Revised parameters: ", err_exp_3_rev_temp_acslx))
}