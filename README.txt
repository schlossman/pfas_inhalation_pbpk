PBPK Model Template Project - pfas_inhalation branch/application

-------------------------
PFOA_inhalation_sim_combinedExpo.R
This is the master script to generate the PFAS PBPK model results shown in the 
U.S. EPA paper (Dzierlenga et al., 2025), "Challenges for extrapolation of risk 
from ingestion to inhalation exposure for per- and polyfluorinated alkyl acids 
and their precursors." Running the script by opening and selecting "Source" in 
RStudio should generate the plots shown in Figures 3, 4 and S1-S3.

-------------------------
Additional documentation including information on how to run these files can be 
found in the Documentation folder.

**The current format for input parameter spreadsheets can be seen in:
model parameters: BLANK_template_parameters_Model.xlsx
exposure parameters: BLANK_template_parameters_Exposure.xlsx


To use the included R functions, GNU MCSim must be installed and the "mod.exe" 
utility must be available in the user's PATH.

-------------------------
run_template_model.R

This R source code file includes functions to run the PBPK model template. 
It requires the R libraries readxl (for importing parameter values) and 
deSolve (for methods of solving ordinary differential equations) and the 
R file RMCSim.R. The included functions are described in 
PBPK_Template_Description_of_Functions.

-------------------------
PBPK_template.model
This MCSim source code file defines the current model template.

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
RMCSim.R
This R source code file defines general functions for compiling, loading, and 
running an ODE model encoded in the GNU MCSim model specification language. 
The mod.exe utility must be available in the user's PATH for the compile_model 
function to work properly.

-------------------------
exposure_functions.R
This R source code file defines functions pulse.inh.exp() which creates a data
frame defining periodic inhalation exposure scenarios and code to define the
data frame for bolus oral dosing events, df_dose.

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
Data
This is a directory containing data files for each of the chemicals and example 
models.  The directory contains sub-folders for each chemical named 
"Digitized_Data_X" where X is the name of the chemical. The sub-folders 
contain .csv and .xlsx files containing data digitized from the source 
publications.

-------------------------
Inputs
This is a directory containing input data files for each chemical. To run the 
current model template there must be two input spreadsheets: one for parameters 
in the model, and one for parameters related to the exposure scenario. Note, 
the top of this file indicates an example file containing the current format 
for the input spreadsheets. Previously implemented chemical models may need to have their input spreadsheets modified in order to be compatible with the current version of the template.