# PFAS Inhalation Model
#rm(list = ls()) # tidy
# Set working directory to the directory containing this file.
script.dir = dirname(sys.frame(1)$ofile)
setwd(script.dir)

old.par <- par(no.readonly = TRUE)
source("run_template_model.R")
source("PFOA_Loccisano_Kemper.R")
# Will use Loccisano et al. (2012) PFOA model parameters except those changed in adj_parms.

# PFOA inhalation study - Hinderliter et al. (2006), DuPont report (2003)
# Load single-exposure data (male rats)
data.loc <- "Data/Data_PFOA_inhalation/"
data.sing.exp <- read.csv(file = paste0(data.loc, "Hinderliter_Dupont_Table6.csv"), header = TRUE, sep = ",")
doses=c(25,10,1)
times = seq(from=0, to=32, by=0.1)
rtol=1e-12
atol=1e-12

# From Longhurst et al. (1992), https://doi.org/10.1016/S0022-5347(17)36981-1,
# the average daily urination rate in 6-month male F344 rats, with an average
# BW of 0.363 kg, was 6.6 mL = 0.0066 L/d. Using the filtrate volume fraction 
# of 0.00084 L/kg, the daily urine production in the model corresponds to 
# r_urine (L/h) = k_ust*V_fil = (k_ustc/BW^0.25)*V_filc*BW = k_ustc*V_filc*BW^0.75. 
# Setting this equal to 0.0066/24 yields: ,
# k_ustc = (0.0066/24)/((BW^0.75)*V_filc) 
#        = (0.0066/24)/((0.363^0.75)*V_filc) = 0.700041 kg^0.25/h

model.sheetname = "MaleRat_2blood_noapprox" 
  # Parameter sheet for two blood compartments, no SS approximation

cases = list(
  case1=list(adj_parms=c(P_B=100000), 
          labysr=c(46,25,12), 
          simtitle="No Published Parameters Changed, Lung & Second Blood Compartments Added"),
  case2=list(adj_parms=c(V_max_reabsc=0.5, 
                      K_m_reabs=60, 
                      k_ustc=0.7, 
                      P_rb=0.04, 
                      k_bilec=0.22, 
                      k_fecesc=0.005, 
                      P_B=100000, 
                      k_unabs=2),
          labysr=c(52,25,12), 
          simtitle="P_RB & Excretion Parameters Adjusted" )
  )

for (case in cases){
  
PFOA.Loccisano.Kemper(case=case, colr=FALSE, model.sheetname=model.sheetname)

if (is.null(case$adj_parms)) {
  case$simtitle = "Single and Multi-Day PFOA Inhalation PK Data & Simulations"
  }
#case$adj_parms = c(case$adj_parms, Q_cardiacc=9.94, VPR=5.18)

exp.sheetname.head = "inhal_oral_6hr_" # Exposure sheet name header for combined inhalation/oral simulation
par(mfrow=c(1,2), mar=c(2.4,2.7,1.75,0.2), oma=c(1.5,0,0,0), mgp=c(0.5,0.5,0))
par(mfrow=c(1,2), mar=c(2.4,2.7,1.75,0.2), oma=c(1.5,0,0,0), mgp=c(1.5,0.5,0))
ltys=c(1,3,2)
pchs=c(1,2,6)
labys=c(18,9.5,3.5)
fils=c()
for (i in 1:length(doses)) {
  out <- PBPK_run(model = template,
                  model.param.filename = "PFOA_inh_template_parameters_Model.xlsx", 
                  model.param.sheetname = model.sheetname, 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = paste0(exp.sheetname.head,doses[i]),
                  rtol=rtol, atol=atol, adj.parms=case$adj_parms)
  time1=out[,"time"]
  q1=out[,"Q_cc"]
  data = subset.data.frame(data.sing.exp, Dose==doses[i])
  if (i==1){
    ymax = max(data$Plasma_Conc_Mean, out[,"C_ven"])
    plot(1,1, type="n", xlab="Time (h)", title=NULL, 
         ylab=expression(paste("Venous blood concentration ( ",mu,"g/mL)")), 
         ylim=c(0,ymax), xlim = c(0,max(time1)))
  }
  points(data$Time, data$Plasma_Conc_Mean, pch=pchs[i])
  lines(time1, out[,"C_ven"],lty=ltys[i])
  if (i==1) {
    text(17,labys[i],expression("25 mg/m"^3))
  } else if (i==2) {
    text(17,labys[i],expression("10 mg/m"^3))
  } else {
    text(17,labys[i],expression("1 mg/m"^3))
  }
  fils=cbind(fils,out[,"C_fil"])
}

# Repeat-dose data & sims
data.rep.exp <- read.csv(file = paste0(data.loc, "Hinderliter_Dupont_Table8.csv"), header = TRUE, sep = ",")
BW.table = PFOA_BW_v_age(BW0=0.29, sex="male", BWdata="NTP_BW_SD.csv", dur=26)
exp.sheetname.head = "inhal_6hr_" # Exposure sheet name header [for combined inhalation/oral simulation]
#par(mar=c(2.7,3,0.75,0.2),mgp=c(1.8,0.7,0))
par(mar=c(2.4,1.7,1.75,0.2))
filr=c()
for (i in 1:length(doses)) {
  out <- PBPK_run(model = template, adj.parms = case$adj_parms,
                  model.param.filename = "PFOA_inh_template_parameters_Model.xlsx", 
                  model.param.sheetname = model.sheetname, 
                  exposure.param.filename = "PFOA_template_parameters_Exposure.xlsx", 
                  exposure.param.sheetname = paste0(exp.sheetname.head,doses[i],"_repeat"),
                  BW.table=BW.table, rtol=rtol, atol=atol, method="bdf_d")
  data = subset.data.frame(data.rep.exp, Dose==doses[i])
  if (i==1){
    #plot(out$time,out$R_oral,type="l")
    #par(mar=c(2.7,3,0.75,0.2),mgp=c(1.8,0.7,0))
    ymax = max(data$Plasma_Conc_Mean, 80)#out$C_ven) #
    plot(1,1, type="n", xlab="Time (h)", ylim=c(0,ymax), xlim=c(0,max(out[,"time"])), 
         ylab=expression(paste("Venous blood concentration (",mu,"g/mL)")), 
         title=NULL)
  }
  points(data$Time, data$Plasma_Conc_Mean, pch=pchs[i])
  lines(out[,"time"], out[,"C_ven"],lty=ltys[i])
  if (i==1) {
    text(250,case$labysr[i],expression("25 mg/m"^3))
  } else if (i==2) {
    text(250,case$labysr[i],expression("10 mg/m"^3))
  } else {
    text(250,case$labysr[i],expression("1 mg/m"^3))
  }
  filr=cbind(filr,out[,"C_fil"])
}
#plot(out$time,out$Q_cc,type="l")
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
title(case$simtitle,line=-1.25)

#Plot filtrate concentrations
par(mfrow=c(1,2), mar=c(2.4,2.7,1.7,0.2), oma=c(0.5,0,0,0), mgp=c(1.3,0.5,0))
plot(time1, fils[,1], mgp=c(1.4,0.5,0), title = NULL, col="red", type = "l", 
     xlab="Time (h)", ylab=expression(paste("Filtrate concentration (",mu,"g/mL)")))
lines(time1, fils[,2],lty=3, col="red") 
lines(time1, fils[,3],lty=2, col="red")

plot(out$time, filr[,1], mgp=c(1.4,0.6,0), title = NULL,
     type = "l", xlab = "Time (hr)", ylab = NA, col="red")
lines(out$time, filr[,2],lty=3, col="red") 
lines(out$time, filr[,3],lty=2, col="red")
par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
title(case$simtitle,line=-1.2)

#plot(time1,q1,type="l")

} # End of 'case in cases' set
