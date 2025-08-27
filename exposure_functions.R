# Exposure functions

pulse.inh.exp <- function(dose, time.exp.starts, length.exp.day, N.days.exp, sim.days){
  # Creates a data frame containing times and changes to the ambient air 
  # concentration (Conc) of chemical, updated in the model using replace method (rep)
  
  # dose is the inhaled concentration for an open chamber (ppm)
  # time.exp.starts is the time the exposure starts each day
  # length.exp.day is the length of time the exposure lasts each day
  # N.days.exp is the number of days per week the exposure occurs, 
  #        assumed to be the first N days of the week
  # sim.days is the number of days in the simulation
  
  # Initialize data frame for dosing events
  df_dose = NULL 
  dose_per_event = c(dose, 0)
  
  # Initialize daily exposure pattern (time on, time off).
  tz_daily = c(time.exp.starts, time.exp.starts+length.exp.day)
  
  # Vector of exposure change times.
  tz = vector("numeric", sim.days * length(tz_daily))
  
  # Build a vector to contain the simulation times at which exposure changes.
  idx = 1
  for (day in seq(1, sim.days, 1)) {
    tz[seq(idx, idx+length(tz_daily)-1, 1)] = (day - 1) * 24 + tz_daily
    idx = idx + length(tz_daily)
  }
  
  # Build a vector containing the indexes of simulation times at which exposure changes.
  keep = N.days.exp*length(tz_daily) # Days with exposure (first N days of week)
  toss = (7-N.days.exp)*length(tz_daily) # Days without exposure (rest of week)
  idx = 1
  toss.idx = c()
  for (week in seq(1, sim.days/7, 1)) {
    toss.idx[seq(idx, idx+toss-1, 1)] = seq(keep+1 + (week-1)*7*length(tz_daily), keep+toss + (week-1)*7*length(tz_daily), 1)
    idx = idx + toss
  }
  
  # Data frame containing bolus dosing events.
  df_dose = data.frame(var=rep("Conc", length(tz)), time=tz,
                       value=rep(dose_per_event), #value=rep(dose_per_event, length(tz)),
                       method=rep("rep", length(tz)))
  # Remove weekend times.
  df_dose = df_dose[-toss.idx,]
  
  # Sorted data frame.
  df_dose = df_dose[order(df_dose$time), ]
  
  return(df_dose)
}

# Code chunk for bolus oral doses (drinking water)
df_dose = NULL
  # if the dose is via drinking water, construct the input
  if(all.exposure.parms$dose.type == "water"){
    # Drinking water dose (mg/d).
    dose = exp.parms[["dose_water"]] * parms[["BW"]]
    
    # Daily water ingestion times: 8am, 10am, 12pm, 2pm, 4pm, and 6pm.
    final.dose = exp.parms[["first.dose_water"]] + exp.parms[["int.doses_water"]]*(exp.parms[["n.doses_water"]]-1)
    tz_daily = seq(from = exp.parms[["first.dose_water"]], to = final.dose, by = exp.parms[["int.doses_water"]])
    
    # Drinking water dose per event (mg).
    if (all.exposure.parms$water.equal == "equal"){
      dose_per_event = dose / length(tz_daily) #equal doses per day
    } else if (all.exposure.parms$water.equal == "unequal") {
      dose_per_event = dose*c(0.233, 0.1, 0.1, 0.1, 0.233, 0.234) #IRIS mouse
    }
    
    
    # Vector of all water dose times.
    tz = vector("numeric", all.exposure.parms$sim.days * length(tz_daily))
    
    # Build a vector to contain the simulation times at which water ingestions
    # occur.
    idx = 1
    for (day in seq(1, all.exposure.parms$sim.days, 1)) {
      tz[seq(idx, idx+length(tz_daily)-1, 1)] = (day - 1) * 24 + tz_daily
      idx = idx + length(tz_daily)
    }
    
    # Data frame containing bolus dosing events.
    df_dose = data.frame(var=rep("A_glumen", length(tz)), time=tz,
                         value=rep(dose_per_event), #value=rep(dose_per_event, length(tz)),
                         method=rep("add", length(tz)))
    
    # Data frame containing updates to total input to system.
    df_in = df_dose
    df_in["var"] = "A_in"
    
    # Sorted data frame.
    df_dose = rbind(df_dose, df_in)
    df_dose = df_dose[order(df_dose$time), ]
    
  } #end if(dose.type == "water")