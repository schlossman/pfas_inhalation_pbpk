PBPK Model Template application: model extrapolation to predict PK from PFOA 
inhalation experiments, with predictions of Kemper et al. oral bolus experiment
for validation, using parameters otherwise from the Loccisano et al. (2012)
PFOA PBPK model.

-------------------------
PFOA_inhalation_sim_combinedExpo.R
This is the *master script* to generate the PFAS PBPK model results shown in 
the U.S. EPA paper (Dzierlenga et al., 2025), "Challenges for extrapolation of 
risk from ingestion to inhalation exposure for per- and polyfluorinated alkyl 
acids and their precursors." Running this script by opening it and selecting 
"Source" in RStudio should generate the plots shown in Figures 3, 4 and S1-S3.

-------------------------
Additional documentation including information on how to run these files can be 
found in the Documentation folder.

**The input parameter spreadsheets can be seen in the "Inputs" folder:
model parameters: PFOA_inh_template_parameters_Model.xlsx
exposure parameters: PFOA_template_parameters_Exposure.xlsx


To use the included R functions, the R packages readxl, MCSimMod and deSolve 
are needed. (Other, dependent packages may install with these.)

-------------------------
PFOA_inhalation_sim_combinedExpo.R

This is a master R script which, when run, will create the MCSimModel PBPK 
Template Model object and source other sccripts and call their functions as 
needed to generate Figures 3, 4 and S1-S3. Additional figures are generated 
showing the concentration of PFOA in the renal filtrate for diagnostic purposes.
Note: some of the figure titles may from those shown in the paper but inform the
viewer whether published parameters are used or (some) are revised to better 
match the data as described in the paper.

-------------------------
PBPK_template.model
This MCSim source code file defines the current model template.

The following files are created by the following commands, in the order shown.
  template <- createModel("PBPK_template")
  template$loadModel()
These are run, respectively, when run_template_model.R is sourced and the 
function PBPK_run() is called for the first time, which happens automatically 
when the master script, PFOA_inhalation_sim_combinedExpo.R, is called.

PBPK_template_model.c
This C source code file defines the current model template.

PBPK_template_model.dll
This DLL file defines the current model template and can be used directly by 
functions in the R package deSolve.

PBPK_template_model.o
This C object code file defines the current model template.

PBPK_template_model_inits.R
This R source code file, created automatically when the .model file is compiled,
defines initialization functions for the current model template.

-------------------------
run_template_model.R
This R source code file includes functions to run the PBPK model template. 
It uses the R libraries readxl (for importing parameter values) and deSolve 
(for methods of solving ordinary differential equations). The defined functions 
are described in PBPK_Template_Description_of_Functions.

-------------------------
optimization_scripts.R
Contains R functions to perform model sensitivity analyses and parameter 
estimation. These are not used for the pfas inhalation analysis conducted in
2025 but the file is retained as a part of the PBPK Model Template project.

-------------------------
PFOA_Loccisano_Kemper.R
Functions to run Loccisano et al. (2012) PFOA PBPK model simulations of Kemper 
et al. oral PK data and plot results.

-------------------------
Data [folder]
This is a directory (with sub-folders) containing data files used. The master
folder contains .csv files with BW data for SD and Wistar rats from NTP and 
an .xlsx file with Wistar rat data. The sub-folders contain .csv and .xlsx 
files containing data digitized from the source publications.

-------------------------
Inputs [folder]
This is a directory containing input data files for each chemical. To run the 
current model template there must be two input spreadsheets: one for parameters 
in the model, and one for parameters related to the exposure scenario. Note, 
previously implemented chemical models (not included) may need to have their 
scripts revised to work with MCSimMod and may need to have their input 
spreadsheets modified in order to be compatible with this version of the PBPK
Model Template.