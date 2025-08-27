# Functions to run simulations for Loccisano et al. (2012) PFOA PBPK model, 
# simulations of Kemper et al. oral PK data, and plot results:
# Paul Schlosser, U.S. EPA, July 2025
# Requires source of "run_template_model.R".

# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

old.par <- par(no.readonly = TRUE)
source("run_template_model.R")

PFOA_BW_v_age <- function(age0=NULL, BW0=NULL, sex="male", BWdata="NTP_BW_SD.csv", 
                          dur=21){
  # If specified, age0 and dur[ation] should be in days, BW0 in kg.
  if (is.null(age0) & is.null(BW0)){
    print("Either age0 or BW0 must be provided as an input.")
    return(NULL)
  }
  d = read.csv(paste("Data/",BWdata,sep=""))
  d[,2:3] = d[,2:3]/1000  # Convert BWs to g
  if (is.numeric(BW0)){
    if (BW0<d[1,sex]|BW0>max(d[,sex])) stop("Input BW0 is outside the range of BWdata being interpolated.")
    age0 = approx(x=d[1:22,sex],y=d$age[1:22],xout=BW0)$y
  }
  ages = c(age0,d$age[(d$age>age0)&(d$age<(age0+dur))],(age0+dur))
  BW = approx(x=d$age,y=d[,sex],xout=ages,rule=2)$y
    # Rule=2 means that final value in d[,sex] is used for times beyond the data file.
  res = list(times = (ages-age0)*24, BW=BW)
  return(res)
}

PFOA.Loccisano.Kemper <- function(img.name = NULL, case=list(), colr=TRUE,
                                  model.sheetname = "MaleRat_2blood_noapprox"){
  # Figure 10: Loccisano 2012, Male rat, 25 mg/kg PFOA
  # Further below, simulations for 5 and 1 mg/kg PFOA 
  
  case$adj_parms = c(case$adj_parms, Q_kic=0.141, Q_rbc=0.676) 
      # Reset kidney/rbc flows to Loccisano 2012 values
  
  # Construct BW table
  BW0 = 0.2098 # From Table 9 of Kemper report for male rats @ 25 mg/kg
  BW.table = PFOA_BW_v_age(BW0=BW0, sex="male", BWdata="NTP_BW_SD.csv", dur=26)
    # Uses reported growth data from NTP for SD rats, but transposed to match
    # the reported initial BW.
  
  out <- PBPK_run(model.param.filename = "PFOA_inh_template_parameters_Model.xlsx", 
                  model.param.sheetname = model.sheetname, 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKemperOral25BW", # dose = 25 mg/kg
                  BW.table = BW.table, adj.parms = case$adj_parms)

  plot.Kemper.Loccisano(out, route="gavage", dose=25.0, sex="Male", 
                        img.name=img.name, colr=colr, detail=case$simtitle)
  print(paste("Maximum mass balance error:", max(abs(out$A_bal))))
  par(mar=c(3,2.4,0.75,0.2), mgp=c(1.5,0.5,0))
  fil=out$C_fil
  
  # Accuracy calculation vs. digitized simulation results from Loccisano
  data.loc <- "Data/Digitized_Data_PFOA/"
  Pdata <- read.csv(file = paste0(data.loc, "Data_Fig9_Kemper_OraldosePlasma.csv"), header = TRUE, sep = ",")
  data.times <- Pdata[,1]
  out.inc.data <- PBPK_run(model.param.filename = "PFOA_inh_template_parameters_Model.xlsx", 
                           model.param.sheetname = model.sheetname, 
                           exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                           exposure.param.sheetname = "MKemperOral25BW", BW.table=BW.table, 
                           data.times=data.times, adj.parms=case$adj_parms)
  out.inc.data <- out.inc.data[-c(1), ] #remove the zero row that was added in PBPK_run()
  perc <- perc.diff(model = out.inc.data$C_ven, data = Pdata[,2])
  # 
  #plot(data.times, perc, type = "p", xlab = "Time (h)", ylab = "% Difference")
  #lines(data.times, data.times*0)
  #perc.scale <- perc.diff.scale(model = out.inc.data$C_ven, data = Pdata[,2], fig.scale = 200)
  #points(data.times, perc.scale, col = "red")

  
  # Construct BW table
  BW0 = 0.1976 # From Table 9 of Kemper report for male rats @ 5 mg/kg
  BW.table = PFOA_BW_v_age(BW0=BW0, sex="male", BWdata="NTP_BW_SD.csv", dur=21)
  
  out <- PBPK_run(model.param.filename = "PFOA_inh_template_parameters_Model.xlsx", 
                  model.param.sheetname = model.sheetname, 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = "MKemperOral5BW",  # 5 mg/kg dose
                  BW.table=BW.table, adj.parms=case$adj_parms)
  
  #plot.Kemper.Loccisano(out, route="gavage", dose=5.0, sex="Male", 
  #                      img.name=img.name, colr=colr)
  
  BW0 = 0.2334 # From Table 9 of Kemper report for male rats @ 1 mg/kg
  BW.table = PFOA_BW_v_age(BW0=BW0, sex="male", BWdata="NTP_BW_SD.csv", dur=21)
  
  out2 <- PBPK_run(model.param.filename = "PFOA_inh_template_parameters_Model.xlsx", 
                   model.param.sheetname = model.sheetname,
                   exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                   exposure.param.sheetname = "MKemperOral1BW",  # 1 mg/kg dose 
                   BW.table = BW.table, adj.parms=case$adj_parms)
  
  plot.Kemper.Loccisano(out, route="gavage", dose=5, sex="Male", detail=case$simtitle, 
                        img.name=img.name, colr=colr, out2=out2, dose2=1.0)
  
  par(old.par)
  plot(out$time.hr, fil, col="red", type="l", ylab = "Concentraton (ug/mL)",
       xlab = "Time (h)")
  lines(out$time.hr, out$C_fil, col="red", lty=3)
  lines(out$time.hr, out2$C_fil, col="red", lty=2)
  title("Renal filtrate concentration predictions")
  title(case$simtitle, line=-1)
}

plot.Kemper.Loccisano <- function(out, route=NULL, dose=NULL, sex=NULL, detail=NULL, 
                                  img.name=NULL, colr=TRUE, dose2=NULL, out2=NULL){
  # Saves as a tif with name img.name if option is included.
  # Note, you must include the ".tif" in the image name 
  
  data.loc <- "Data/Digitized_Data_PFOA/"
  d <- read.csv(file=paste0(data.loc,"Kemper_6302380.csv"), header=TRUE, sep=",")
  plasma.data <- d[(d$dose==dose)&(d$route==route)&(d$sex==sex),c("time","conc_mean")]
  
  # Set colors from the viridis palette. The first vector of color-identification
  # strings was generated with the function "viridis" from the package "viridis".
  
  if (colr) {
    vcol = c("#440154FF", "#2A788EFF", "#7AD151FF")
  } else { # For grey-scale plots
    vcol = c("lightgrey","black","darkgrey")
  }
  
  # Set colors and line-types for results from Loccisano and PBPK Template:
  paper.col <- vcol[2]
  inc.templ.col <- vcol[1]
  cor.templ.col <- vcol[3]
  
  paper.lty <- "dashed"
  inc.templ.lty <- "dotted"
  cor.templ.lty <- "solid"
  
  # Plotting labels:
  ylabel = "PFOA concentration (ug/mL)"
  xlabel = "Time (h)"
  
  if (!is.null(img.name)){ # Create tiff if file name is not null.
    tiff(img.name, res=300, height=4, width=7, units="in")
  }
  # Set plot frame parameters:
  opar <- par(no.readonly = TRUE)
  par(mfrow=c(1,2), mar=c(3,2.4,1.75,0.2), oma=c(1.5,0.1,0.1,0.1), mgp=c(1.5,0.5,0))
  if (dose==25){ 
    par(mar=c(0.2,2.4,1.75,0.2), oma=c(0.1,0.1,0.1,0.1), xaxt="n")
    xlab=NA
  } else {
    par(mar=c(2.7,2.4,0,0.2))
  }
  xmax <- max(out$time)*1.01
  if (out$C_ven[4]==0) {blood=out$C_bl} else {blood=out$C_ven}
  y_min <- min(blood[blood > 0], plasma.data[plasma.data[,2]>0,2], na.rm = TRUE)
  y_max <- max(blood, plasma.data[,2], na.rm = TRUE)
  Pylim <- c(y_min, y_max)  # y-axis limits for plotting
  
  # Create plot for plasma concentration vs time.
  plot(1,1, type = "n", xlab = xlabel, xlim = c(0,xmax), # Create (sub)plot
       ylab = ylabel, ylim = Pylim, log = "y")
  points(plasma.data[,1],plasma.data[,2], pch = 19, col = paper.col) # Plot Kemper data
  lines(out$time.hr, blood, lty = cor.templ.lty, col = cor.templ.col, lwd = 3) # Plot simulation results
  title(paste("Plasma,",dose,"mg/kg"), line=-1)
  
  if (dose==25) { # Create 2nd plot for excretion data
    # Get (published) Loccisano model simulation results & excretion data (male rats, 25 mg/kg dose)
    pmodel <- read.csv(file = paste0(data.loc, "Data_Fig9_Kemper_OraldosePlasma.csv"), header = TRUE, sep = ",")[,1:2]
    excr.data <- read.csv(file = paste0(data.loc, "Data_Fig9_Kemper_OraldoseExcretion.csv"), header = TRUE, sep = ",")
    lines(pmodel[,1], pmodel[,2], lty = paper.lty, col = paper.col, lwd = 2)
    
    # Convert PBPK Template excretion results to percents.
    urine.percent <- 100*out$A_urine/tail(out$A_in,1)
    feces.percent <- 100*out$A_fecal/tail(out$A_in,1)
    exhale.percent <-  100*out$A_exh/tail(out$A_in,1)
    
    y_min <- min(excr.data[,2], excr.data[,6], na.rm = TRUE)
    y_max <- max(excr.data[,2], excr.data[,4], excr.data[,6], urine.percent, feces.percent, na.rm = TRUE)
    Eylim <- c(y_min,y_max) # y-axis limits for excretion plotting
    
    ylabel <- "Cumulative Percent of Dose"
    
    # Create plot for excretion data and exhaled amount vs time:
    plot(1,1, type = "n", xlab = xlabel, xlim = c(0,xmax),
         ylab = ylabel, ylim = Eylim, log = "y")
    # Urine results:
    lines(out$time.hr, urine.percent, lty = cor.templ.lty, col = cor.templ.col, lwd = 3)
    points(excr.data[,3], excr.data[,4], pch = 19, col = paper.col)
    lines(excr.data[,1], excr.data[,2], lty = paper.lty, col = paper.col, lwd = 2)
    text(500,0.9*excr.data[excr.data[,1]>=500,2][1],"Urine",font=2)
    
    # Fecal results:
    lines(out$time.hr, feces.percent, lty = cor.templ.lty, col = cor.templ.col, lwd = 3)
    points(excr.data[,7], excr.data[,8], pch = 19, col = paper.col)
    lines(excr.data[,5], excr.data[,6], lty = paper.lty, col = paper.col, lwd = 2)
    tht = max(excr.data[excr.data[,5]>=500,6][1],feces.percent[out$time.hr==500])
    text(500,1.15*tht,"Feces",font=2)
    
    # Exhalation results
    lines(out$time.hr, exhale.percent, lty = "dotted", col = cor.templ.col, lwd = 3)
    text(500,1.15*exhale.percent[out$time.hr==500],"Exhaled Air",font=2)
    
    title(paste("Excretion,",dose,"mg/kg"), line=-1)
  } 
  
  if (!is.null(dose2)){ # Assumes dose is not 25, so there's still room for a
                        # 2nd sub-plot. Assumes sex == "male".
                        # Plot plasma results for dose2. dose2 should be 5 or 1.
    # Corresponding plasma data:
    plasma.data <- d[(d$dose==dose2)&(d$route==route)&(d$sex==sex),c("time","conc_mean")]
    xmax <- max(out2$time)*1.01 
    # Get sim data depending on whether 1 or 2 blood compartments used
    if (out2$C_ven[4]==0) {blood=out2$C_bl} else {blood=out2$C_ven}
    y_min <- min(blood[blood > 0], plasma.data[plasma.data[,2]>0,2], na.rm = TRUE)
    y_max <- max(blood, plasma.data[,2], na.rm = TRUE)
    Pylim <- c(y_min, y_max) # y-axis limits for plotting
    
    plot(1,2, type = "n", xlab = xlabel, xlim = c(0,xmax),
         ylab = ylabel, ylim = Pylim, log = "y")
    points(plasma.data[,1], plasma.data[,2], pch = 19, col = paper.col)
    
    lines(out$time.hr, blood, lty = cor.templ.lty, col = cor.templ.col, lwd = 3)
    
    title(paste("Plasma,",dose2,"mg/kg"),line=-1)
  }
  
  par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  
  if (dose==25) {
    title(detail,line=-1.25)
  } else {
    legend("bottom", inset=c(0, 0), xpd=TRUE, bty="n", lwd=c(1,2,3), horiz=TRUE,  
           legend = c("Published Data  ", "Published Model ", "Template Predictions"),
           col = c(paper.col, paper.col, cor.templ.col), 
           lty = c(NA, paper.lty, cor.templ.lty), pch = c(19,NA,NA))}
  if (!is.null(img.name)){ dev.off() }
  par(opar)
}



