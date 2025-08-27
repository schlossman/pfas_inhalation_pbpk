/* PBPK_template_model.c for R deSolve package
   ___________________________________________________

   Model File:  PBPK_template.model

   Date:  Wed May 28 14:47:36 2025

   Created by:  "C:/Users/pschloss/ONEDRI~1/MCSIM-~1.5/mod/mod.exe v6.1.0"
    -- a model preprocessor by Don Maszle
   ___________________________________________________

   Copyright (c) 1993-2019 Free Software Foundation, Inc.

   Model calculations for compartmental model:

   54 States:
     A_bl = 0.0,
     A_ven = 0.0,
     A_art = 0.0,
     A_lu = 0.0,
     A_gi = 0.0,
     A_li = 0.0,
     A_fst = 0.0,
     A_ki = 0.0,
     A_fil = 0.0,
     A_ust = 0.0,
     A_tc1 = 0.0,
     A_tc2 = 0.0,
     A_tc3 = 0.0,
     A_tc4 = 0.0,
     A_tc5 = 0.0,
     A_om = 0.0,
     A_rb = 0.0,
     A_glumen = 0.0,
     A_lib = 0.0,
     R_IV = 0.0,
     R_oral = 0.0,
     A_in = 0.0,
     A_loss = 0.0,
     A_exh = 0.0,
     A_urine = 0.0,
     A_fecal = 0.0,
     A_ch = 0.0,
     A_met_sat = 0.0,
     A_met_sat_li = 0.0,
     A_met_sat_om = 0.0,
     A_met_sat_lu = 0.0,
     A_met_1st = 0.0,
     A_met_1st_li = 0.0,
     A_met_1st_lu = 0.0,
     A_met_1st_om = 0.0,
     V_max_om_t = 0.0,
     V_max_li_t = 0.0,
     V_max_lu_t = 0.0,
     AUC_bl = 0.0,
     AUC_ven = 0.0,
     AUC_art = 0.0,
     AUC_lu = 0.0,
     AUC_gi = 0.0,
     AUC_li = 0.0,
     AUC_ki = 0.0,
     AUC_tc1 = 0.0,
     AUC_tc2 = 0.0,
     AUC_tc3 = 0.0,
     AUC_tc4 = 0.0,
     AUC_tc5 = 0.0,
     AUC_om = 0.0,
     AUC_rb = 0.0,
     Conc = 0.0,
     Q_cc = 0.0,

   105 Outputs:
    "A_bal",
    "A_body",
    "A_out",
    "C_li",
    "C_lu",
    "C_gi",
    "C_ki",
    "C_fil",
    "C_tc1",
    "C_tc2",
    "C_tc3",
    "C_tc4",
    "C_tc5",
    "C_om",
    "C_rb",
    "C_bl",
    "C_ch",
    "C_chppm",
    "C_inh",
    "C_ven",
    "C_pulv",
    "C_art_in",
    "C_art_comp",
    "C_art",
    "BW_out",
    "Free",
    "CarD",
    "Q_bal",
    "Q_cardiac",
    "Q_li",
    "Q_ki",
    "Q_gi",
    "Q_tc1",
    "Q_tc2",
    "Q_tc3",
    "Q_tc4",
    "Q_tc5",
    "Q_om",
    "Q_rb",
    "k_GFR",
    "Q_p",
    "V_max_reabs",
    "k_bile",
    "k_urine",
    "k_feces",
    "k_ust",
    "k_bloss",
    "V_max_lu",
    "k_met_lu",
    "V_max_li",
    "k_met_li",
    "V_max_om",
    "k_met_om",
    "R_met_sat_li",
    "R_met_1st_li",
    "R_met_sat_om",
    "R_met_1st_om",
    "R_met_sat_lu",
    "R_met_1st_lu",
    "R_artfil",
    "R_libind",
    "R_reabs",
    "R_blloss",
    "R_glgi",
    "R_glli",
    "R_glfst",
    "R_gifst",
    "R_gili",
    "R_biliary",
    "R_liunbind",
    "R_feces",
    "R_filust",
    "R_urine",
    "R_exh",
    "R_inh",
    "V_bl",
    "V_ven",
    "V_art",
    "V_lu",
    "V_ki",
    "V_fil",
    "V_li",
    "V_gi",
    "V_tc1",
    "V_tc2",
    "V_tc3",
    "V_tc4",
    "V_tc5",
    "V_om",
    "V_rb",
    "V_total_bl_value",
    "V_lu_value",
    "V_ki_value",
    "V_fil_value",
    "V_li_value",
    "V_gi_value",
    "V_tc1_value",
    "V_tc2_value",
    "V_tc3_value",
    "V_tc4_value",
    "V_tc5_value",
    "V_om_value",
    "V_rb_value",
    "Q_relative",
    "total_volume_frac",

   2 Inputs:
     BW_in (forcing function)
     Free_in (forcing function)

   101 Parameters:
     MOLWT = 0.0,
     RTemp = 24450.0,
     BW = 0.0,
     F_free = 1.0,
     oral_dose_init = 0.0,
     iv_dose = 0.0,
     single_blood = 0.0,
     venous_ss = 1.0,
     arterial_ss = 1.0,
     exist_lung = 1.0,
     GE_ss = 1.0,
     Conc_init = 0.0,
     NCH = 0.0,
     VCHC = 0.0,
     KL = 0.0,
     F_inh = 1.0,
     AS_co = 0.75,
     AS_met = 0.75,
     AS_cl_met = -0.25,
     AS_cl = -0.25,
     P_lu = 1.0,
     P_li = 1.0,
     P_ki = 1.0,
     P_gi = 1.0,
     P_tc1 = 1.0,
     P_tc2 = 1.0,
     P_tc3 = 1.0,
     P_tc4 = 1.0,
     P_tc5 = 1.0,
     P_om = 1.0,
     P_rb = 1.0,
     V_max_reabsc = 0.0,
     K_m_reabs = 1.0,
     k_urinec = 0.0,
     k_fst = 0.0,
     k_absgi = 0.0,
     k_absli = 0.0,
     k_absli2 = 0.0,
     F_unabs = 0.0,
     k_unabs = 0.0,
     k_ustc = 0.0,
     k_ven_ustc = 0.0,
     k_bilec = 0.0,
     k_fecesc = 0.0,
     k_loss = 0.0,
     VPR = 0.0,
     P_B = 1.0,
     V_max_bind_li = 0.0,
     K_m_bind_li = 1.0,
     k_off_li = 0.0,
     k_blossc = 0.0,
     k_GFRc = 0.0,
     k_enz_resyn = 0.0,
     k_enz_loss = 0.0,
     V_max_luc = 0.0,
     K_m_lu = 1.0,
     k_met_luc = 0.0,
     V_max_lic = 0.0,
     K_m_li = 1.0,
     k_met_lic = 0.0,
     V_max_omc = 0.0,
     K_m_om = 1.0,
     k_met_omc = 0.0,
     Q_cardiacc = 0.0,
     Q_lic = 0.0,
     Q_kic = 0.0,
     Q_gic = 0.0,
     Q_tc1c = 0.0,
     Q_tc2c = 0.0,
     Q_tc3c = 0.0,
     Q_tc4c = 0.0,
     Q_tc5c = 0.0,
     Q_omc = 0.0,
     Q_rbc = 0.0,
     V_blc = 1.0,
     V_venc = 1.0,
     V_artc = 1.0,
     V_luc = 1.0,
     V_kic = 1.0,
     V_filc = 1.0,
     V_lic = 1.0,
     V_gic = 1.0,
     V_tc1c = 1.0,
     V_tc2c = 1.0,
     V_tc3c = 1.0,
     V_tc4c = 1.0,
     V_tc5c = 1.0,
     V_omc = 1.0,
     V_rbc = 1.0,
     R_0bgli = 0.0,
     CDSW = 0.0,
     KIV = 1.0,
     NV = 0.0,
     GASD = 0.0,
     V_ch = 0.0,
     A_bl0 = 0.0,
     A_ven0 = 0.0,
     Conc_ambient = 0.0,
     V_max_om_0 = 0.0,
     V_max_li_0 = 0.0,
     V_max_lu_0 = 0.0,
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Model variables: States */
#define ID_A_bl 0x00000
#define ID_A_ven 0x00001
#define ID_A_art 0x00002
#define ID_A_lu 0x00003
#define ID_A_gi 0x00004
#define ID_A_li 0x00005
#define ID_A_fst 0x00006
#define ID_A_ki 0x00007
#define ID_A_fil 0x00008
#define ID_A_ust 0x00009
#define ID_A_tc1 0x0000a
#define ID_A_tc2 0x0000b
#define ID_A_tc3 0x0000c
#define ID_A_tc4 0x0000d
#define ID_A_tc5 0x0000e
#define ID_A_om 0x0000f
#define ID_A_rb 0x00010
#define ID_A_glumen 0x00011
#define ID_A_lib 0x00012
#define ID_R_IV 0x00013
#define ID_R_oral 0x00014
#define ID_A_in 0x00015
#define ID_A_loss 0x00016
#define ID_A_exh 0x00017
#define ID_A_urine 0x00018
#define ID_A_fecal 0x00019
#define ID_A_ch 0x0001a
#define ID_A_met_sat 0x0001b
#define ID_A_met_sat_li 0x0001c
#define ID_A_met_sat_om 0x0001d
#define ID_A_met_sat_lu 0x0001e
#define ID_A_met_1st 0x0001f
#define ID_A_met_1st_li 0x00020
#define ID_A_met_1st_lu 0x00021
#define ID_A_met_1st_om 0x00022
#define ID_V_max_om_t 0x00023
#define ID_V_max_li_t 0x00024
#define ID_V_max_lu_t 0x00025
#define ID_AUC_bl 0x00026
#define ID_AUC_ven 0x00027
#define ID_AUC_art 0x00028
#define ID_AUC_lu 0x00029
#define ID_AUC_gi 0x0002a
#define ID_AUC_li 0x0002b
#define ID_AUC_ki 0x0002c
#define ID_AUC_tc1 0x0002d
#define ID_AUC_tc2 0x0002e
#define ID_AUC_tc3 0x0002f
#define ID_AUC_tc4 0x00030
#define ID_AUC_tc5 0x00031
#define ID_AUC_om 0x00032
#define ID_AUC_rb 0x00033
#define ID_Conc 0x00034
#define ID_Q_cc 0x00035

/* Model variables: Outputs */
#define ID_A_bal 0x00000
#define ID_A_body 0x00001
#define ID_A_out 0x00002
#define ID_C_li 0x00003
#define ID_C_lu 0x00004
#define ID_C_gi 0x00005
#define ID_C_ki 0x00006
#define ID_C_fil 0x00007
#define ID_C_tc1 0x00008
#define ID_C_tc2 0x00009
#define ID_C_tc3 0x0000a
#define ID_C_tc4 0x0000b
#define ID_C_tc5 0x0000c
#define ID_C_om 0x0000d
#define ID_C_rb 0x0000e
#define ID_C_bl 0x0000f
#define ID_C_ch 0x00010
#define ID_C_chppm 0x00011
#define ID_C_inh 0x00012
#define ID_C_ven 0x00013
#define ID_C_pulv 0x00014
#define ID_C_art_in 0x00015
#define ID_C_art_comp 0x00016
#define ID_C_art 0x00017
#define ID_BW_out 0x00018
#define ID_Free 0x00019
#define ID_CarD 0x0001a
#define ID_Q_bal 0x0001b
#define ID_Q_cardiac 0x0001c
#define ID_Q_li 0x0001d
#define ID_Q_ki 0x0001e
#define ID_Q_gi 0x0001f
#define ID_Q_tc1 0x00020
#define ID_Q_tc2 0x00021
#define ID_Q_tc3 0x00022
#define ID_Q_tc4 0x00023
#define ID_Q_tc5 0x00024
#define ID_Q_om 0x00025
#define ID_Q_rb 0x00026
#define ID_k_GFR 0x00027
#define ID_Q_p 0x00028
#define ID_V_max_reabs 0x00029
#define ID_k_bile 0x0002a
#define ID_k_urine 0x0002b
#define ID_k_feces 0x0002c
#define ID_k_ust 0x0002d
#define ID_k_bloss 0x0002e
#define ID_V_max_lu 0x0002f
#define ID_k_met_lu 0x00030
#define ID_V_max_li 0x00031
#define ID_k_met_li 0x00032
#define ID_V_max_om 0x00033
#define ID_k_met_om 0x00034
#define ID_R_met_sat_li 0x00035
#define ID_R_met_1st_li 0x00036
#define ID_R_met_sat_om 0x00037
#define ID_R_met_1st_om 0x00038
#define ID_R_met_sat_lu 0x00039
#define ID_R_met_1st_lu 0x0003a
#define ID_R_artfil 0x0003b
#define ID_R_libind 0x0003c
#define ID_R_reabs 0x0003d
#define ID_R_blloss 0x0003e
#define ID_R_glgi 0x0003f
#define ID_R_glli 0x00040
#define ID_R_glfst 0x00041
#define ID_R_gifst 0x00042
#define ID_R_gili 0x00043
#define ID_R_biliary 0x00044
#define ID_R_liunbind 0x00045
#define ID_R_feces 0x00046
#define ID_R_filust 0x00047
#define ID_R_urine 0x00048
#define ID_R_exh 0x00049
#define ID_R_inh 0x0004a
#define ID_V_bl 0x0004b
#define ID_V_ven 0x0004c
#define ID_V_art 0x0004d
#define ID_V_lu 0x0004e
#define ID_V_ki 0x0004f
#define ID_V_fil 0x00050
#define ID_V_li 0x00051
#define ID_V_gi 0x00052
#define ID_V_tc1 0x00053
#define ID_V_tc2 0x00054
#define ID_V_tc3 0x00055
#define ID_V_tc4 0x00056
#define ID_V_tc5 0x00057
#define ID_V_om 0x00058
#define ID_V_rb 0x00059
#define ID_V_total_bl_value 0x0005a
#define ID_V_lu_value 0x0005b
#define ID_V_ki_value 0x0005c
#define ID_V_fil_value 0x0005d
#define ID_V_li_value 0x0005e
#define ID_V_gi_value 0x0005f
#define ID_V_tc1_value 0x00060
#define ID_V_tc2_value 0x00061
#define ID_V_tc3_value 0x00062
#define ID_V_tc4_value 0x00063
#define ID_V_tc5_value 0x00064
#define ID_V_om_value 0x00065
#define ID_V_rb_value 0x00066
#define ID_Q_relative 0x00067
#define ID_total_volume_frac 0x00068

/* Parameters */
static double parms[101];

#define MOLWT parms[0]
#define RTemp parms[1]
#define BW parms[2]
#define F_free parms[3]
#define oral_dose_init parms[4]
#define iv_dose parms[5]
#define single_blood parms[6]
#define venous_ss parms[7]
#define arterial_ss parms[8]
#define exist_lung parms[9]
#define GE_ss parms[10]
#define Conc_init parms[11]
#define NCH parms[12]
#define VCHC parms[13]
#define KL parms[14]
#define F_inh parms[15]
#define AS_co parms[16]
#define AS_met parms[17]
#define AS_cl_met parms[18]
#define AS_cl parms[19]
#define P_lu parms[20]
#define P_li parms[21]
#define P_ki parms[22]
#define P_gi parms[23]
#define P_tc1 parms[24]
#define P_tc2 parms[25]
#define P_tc3 parms[26]
#define P_tc4 parms[27]
#define P_tc5 parms[28]
#define P_om parms[29]
#define P_rb parms[30]
#define V_max_reabsc parms[31]
#define K_m_reabs parms[32]
#define k_urinec parms[33]
#define k_fst parms[34]
#define k_absgi parms[35]
#define k_absli parms[36]
#define k_absli2 parms[37]
#define F_unabs parms[38]
#define k_unabs parms[39]
#define k_ustc parms[40]
#define k_ven_ustc parms[41]
#define k_bilec parms[42]
#define k_fecesc parms[43]
#define k_loss parms[44]
#define VPR parms[45]
#define P_B parms[46]
#define V_max_bind_li parms[47]
#define K_m_bind_li parms[48]
#define k_off_li parms[49]
#define k_blossc parms[50]
#define k_GFRc parms[51]
#define k_enz_resyn parms[52]
#define k_enz_loss parms[53]
#define V_max_luc parms[54]
#define K_m_lu parms[55]
#define k_met_luc parms[56]
#define V_max_lic parms[57]
#define K_m_li parms[58]
#define k_met_lic parms[59]
#define V_max_omc parms[60]
#define K_m_om parms[61]
#define k_met_omc parms[62]
#define Q_cardiacc parms[63]
#define Q_lic parms[64]
#define Q_kic parms[65]
#define Q_gic parms[66]
#define Q_tc1c parms[67]
#define Q_tc2c parms[68]
#define Q_tc3c parms[69]
#define Q_tc4c parms[70]
#define Q_tc5c parms[71]
#define Q_omc parms[72]
#define Q_rbc parms[73]
#define V_blc parms[74]
#define V_venc parms[75]
#define V_artc parms[76]
#define V_luc parms[77]
#define V_kic parms[78]
#define V_filc parms[79]
#define V_lic parms[80]
#define V_gic parms[81]
#define V_tc1c parms[82]
#define V_tc2c parms[83]
#define V_tc3c parms[84]
#define V_tc4c parms[85]
#define V_tc5c parms[86]
#define V_omc parms[87]
#define V_rbc parms[88]
#define R_0bgli parms[89]
#define CDSW parms[90]
#define KIV parms[91]
#define NV parms[92]
#define GASD parms[93]
#define V_ch parms[94]
#define A_bl0 parms[95]
#define A_ven0 parms[96]
#define Conc_ambient parms[97]
#define V_max_om_0 parms[98]
#define V_max_li_0 parms[99]
#define V_max_lu_0 parms[100]

/* Forcing (Input) functions */
static double forc[2];

#define BW_in forc[0]
#define Free_in forc[1]

/* Function definitions for delay differential equations */

int Nout=1;
int nr[1]={0};
double ytau[1] = {0.0};

static double yini[54] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; /*Array of initial state variables*/

void lagvalue(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  return fun(T, nr, N, ytau);
}

double CalcDelay(int hvar, double dTime, double delay) {
  double T = dTime-delay;
  if (dTime > delay){
    nr[0] = hvar;
    lagvalue( T, nr, Nout, ytau );
}
  else{
    ytau[0] = yini[hvar];
}
  return(ytau[0]);
}

/*----- Initializers */
void initmod (void (* odeparms)(int *, double *))
{
  int N=101;
  odeparms(&N, parms);
}

void initforc (void (* odeforcs)(int *, double *))
{
  int N=2;
  odeforcs(&N, forc);
}


/* Calling R code will ensure that input y has same
   dimension as yini */
void initState (double *y)
{
  int i;

  for (i = 0; i < sizeof(yini) / sizeof(yini[0]); i++)
  {
    yini[i] = y[i];
  }
}

void getParms (double *inParms, double *out, int *nout) {
/*----- Model scaling */

  int i;

  for (i = 0; i < *nout; i++) {
    parms[i] = inParms[i];
  }


  GASD = MOLWT / RTemp ;

  V_ch = ( VCHC ? VCHC - ( NCH * BW ) : 0 ) ;

  A_bl0 = ( single_blood ? iv_dose * BW : 0.0 ) ;
  A_ven0 = ( single_blood ? 0.0 : iv_dose * BW ) ;

  Conc_ambient = ( V_ch ? 0.0 : Conc_init ) ;

  V_max_om_0 = V_max_omc * pow ( BW , AS_met ) ;
  V_max_li_0 = V_max_lic * pow ( BW , AS_met ) ;
  V_max_lu_0 = V_max_luc * pow ( BW , AS_met ) ;

  for (i = 0; i < *nout; i++) {
    out[i] = parms[i];
  }
  }
/*----- Dynamics section */

void derivs (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{
  /* local */ double k_ven_ust;
  /* local */ double R_ven_in;
  /* local */ double Ven_loss_c;
  /* local */ double C_ven_ss_lung;
  /* local */ double C_ven_ss;
  /* local */ double C_pulv_ss;
  /* local */ double C_art_ss_lung;
  /* local */ double C_art_ss;
  /* local */ double R_venust;
  /* local */ double R_enz_loss_lu;
  /* local */ double R_enz_loss_li;
  /* local */ double R_enz_loss_om;
  /* local */ double R_ch_loss;
  /* local */ double R_Qgi_in;
  /* local */ double R_Qgi_out;
  /* local */ double R_Qli_in;
  /* local */ double R_Qli_out;
  /* local */ double R_Qki_in;
  /* local */ double R_Qki_out;
  /* local */ double R_Qtc1_in;
  /* local */ double R_Qtc1_out;
  /* local */ double R_Qtc2_in;
  /* local */ double R_Qtc2_out;
  /* local */ double R_Qtc3_in;
  /* local */ double R_Qtc3_out;
  /* local */ double R_Qtc4_in;
  /* local */ double R_Qtc4_out;
  /* local */ double R_Qtc5_in;
  /* local */ double R_Qtc5_out;
  /* local */ double R_Qom_in;
  /* local */ double R_Qom_out;
  /* local */ double R_Qrb_in;
  /* local */ double R_Qrb_out;
  /* local */ double Q_body;

  yout[ID_BW_out] = BW_in ;

  yout[ID_Free] = Free_in ;

  yout[ID_V_max_reabs] = V_max_reabsc * pow ( yout[ID_BW_out] , AS_met ) ;

  yout[ID_V_max_lu] = V_max_luc * pow ( yout[ID_BW_out] , AS_met ) ;

  yout[ID_V_max_li] = V_max_lic * pow ( yout[ID_BW_out] , AS_met ) ;

  yout[ID_V_max_om] = V_max_omc * pow ( yout[ID_BW_out] , AS_met ) ;

  yout[ID_k_met_lu] = k_met_luc * pow ( yout[ID_BW_out] , AS_cl_met ) ;

  yout[ID_k_met_li] = k_met_lic * pow ( yout[ID_BW_out] , AS_cl_met ) ;

  yout[ID_k_met_om] = k_met_omc * pow ( yout[ID_BW_out] , AS_cl_met ) ;

  yout[ID_k_bile] = k_bilec * pow ( yout[ID_BW_out] , AS_cl ) ;

  yout[ID_k_feces] = k_fecesc * pow ( yout[ID_BW_out] , AS_cl ) ;

  yout[ID_k_ust] = k_ustc * pow ( yout[ID_BW_out] , AS_cl ) ;

  k_ven_ust = k_ven_ustc * pow ( yout[ID_BW_out] , AS_cl ) ;

  yout[ID_k_urine] = k_urinec * pow ( yout[ID_BW_out] , AS_cl ) ;

  yout[ID_k_bloss] = k_blossc * pow ( yout[ID_BW_out] , AS_cl ) ;

  yout[ID_V_bl] = V_blc * yout[ID_BW_out] ;

  yout[ID_V_ven] = V_venc * yout[ID_BW_out] ;

  yout[ID_V_art] = V_artc * yout[ID_BW_out] ;

  yout[ID_V_lu] = V_luc * yout[ID_BW_out] ;

  yout[ID_V_ki] = V_kic * yout[ID_BW_out] ;

  yout[ID_V_fil] = V_filc * yout[ID_BW_out] ;

  yout[ID_V_li] = V_lic * yout[ID_BW_out] ;

  yout[ID_V_gi] = V_gic * yout[ID_BW_out] ;

  yout[ID_V_tc1] = V_tc1c * yout[ID_BW_out] ;

  yout[ID_V_tc2] = V_tc2c * yout[ID_BW_out] ;

  yout[ID_V_tc3] = V_tc3c * yout[ID_BW_out] ;

  yout[ID_V_tc4] = V_tc4c * yout[ID_BW_out] ;

  yout[ID_V_tc5] = V_tc5c * yout[ID_BW_out] ;

  yout[ID_V_om] = V_omc * yout[ID_BW_out] ;

  yout[ID_V_rb] = V_rbc * yout[ID_BW_out] ;

  yout[ID_C_bl] = y[ID_A_bl] / yout[ID_V_bl] ;

  yout[ID_C_lu] = y[ID_A_lu] / yout[ID_V_lu] ;

  yout[ID_C_gi] = y[ID_A_gi] / yout[ID_V_gi] ;

  yout[ID_C_li] = y[ID_A_li] / yout[ID_V_li] ;

  yout[ID_C_ki] = y[ID_A_ki] / yout[ID_V_ki] ;

  yout[ID_C_fil] = y[ID_A_fil] / yout[ID_V_fil] ;

  yout[ID_C_tc1] = y[ID_A_tc1] / yout[ID_V_tc1] ;

  yout[ID_C_tc2] = y[ID_A_tc2] / yout[ID_V_tc2] ;

  yout[ID_C_tc3] = y[ID_A_tc3] / yout[ID_V_tc3] ;

  yout[ID_C_tc4] = y[ID_A_tc4] / yout[ID_V_tc4] ;

  yout[ID_C_tc5] = y[ID_A_tc5] / yout[ID_V_tc5] ;

  yout[ID_C_om] = y[ID_A_om] / yout[ID_V_om] ;

  yout[ID_C_rb] = y[ID_A_rb] / yout[ID_V_rb] ;

  yout[ID_CarD] = 1.0 + CDSW * pow ( ( ( yout[ID_C_rb] / P_rb ) / KIV ) , NV ) ;

  yout[ID_Q_cardiac] = y[ID_Q_cc] * pow ( yout[ID_BW_out] , AS_co ) / yout[ID_CarD] ;

  yout[ID_Q_li] = Q_lic * yout[ID_Q_cardiac] ;

  yout[ID_Q_ki] = Q_kic * yout[ID_Q_cardiac] ;

  yout[ID_Q_gi] = Q_gic * yout[ID_Q_cardiac] ;

  yout[ID_Q_tc1] = Q_tc1c * yout[ID_Q_cardiac] ;

  yout[ID_Q_tc2] = Q_tc2c * yout[ID_Q_cardiac] ;

  yout[ID_Q_tc3] = Q_tc3c * yout[ID_Q_cardiac] ;

  yout[ID_Q_tc4] = Q_tc4c * yout[ID_Q_cardiac] ;

  yout[ID_Q_tc5] = Q_tc5c * yout[ID_Q_cardiac] ;

  yout[ID_Q_om] = Q_omc * yout[ID_Q_cardiac] ;

  yout[ID_Q_rb] = Q_rbc * yout[ID_Q_cardiac] ;

  yout[ID_k_GFR] = k_GFRc * yout[ID_Q_cardiac] ;

  yout[ID_Q_p] = VPR * yout[ID_Q_cardiac] ;

  yout[ID_C_ch] = ( V_ch ? y[ID_A_ch] / V_ch : 0 ) ;

  yout[ID_C_chppm] = yout[ID_C_ch] / GASD ;

  yout[ID_C_inh] = ( V_ch ? yout[ID_C_ch] : y[ID_Conc] * GASD ) ;

  R_ven_in = yout[ID_Free] * ( ( yout[ID_Q_li] + yout[ID_Q_gi] ) * yout[ID_C_li] / P_li + yout[ID_Q_ki] * yout[ID_C_ki] / P_ki + yout[ID_Q_tc1] * yout[ID_C_tc1] / P_tc1 + yout[ID_Q_tc2] * yout[ID_C_tc2] / P_tc2 + yout[ID_Q_tc3] * yout[ID_C_tc3] / P_tc3 + yout[ID_Q_tc4] * yout[ID_C_tc4] / P_tc4 + yout[ID_Q_tc5] * yout[ID_C_tc5] / P_tc5 + yout[ID_Q_om] * yout[ID_C_om] / P_om + yout[ID_Q_rb] * yout[ID_C_rb] / P_rb ) + y[ID_R_IV] ;

  Ven_loss_c = k_ven_ust * yout[ID_V_ven] + k_loss ;

  C_ven_ss_lung = ( GE_ss ? ( ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * R_ven_in + ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_C_lu] ) / P_lu ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) : ( ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) * R_ven_in + ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_C_lu] ) / P_lu ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) ) ;

  C_ven_ss = ( exist_lung ? C_ven_ss_lung : ( ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * R_ven_in + yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * ( yout[ID_Q_p] * F_inh * yout[ID_C_inh] ) ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * pow ( yout[ID_Q_cardiac] * yout[ID_Free] , 2 ) - ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) ) ;

  yout[ID_C_ven] = ( venous_ss ? C_ven_ss : y[ID_A_ven] / yout[ID_V_ven] ) ;

  C_pulv_ss = ( exist_lung ? ( ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) * ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_p] * F_inh * yout[ID_C_inh] ) + ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_C_lu] ) / P_lu + ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) * yout[ID_Q_cardiac] * yout[ID_Free] * R_ven_in - ( yout[ID_Q_p] * F_inh * yout[ID_C_inh] ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) : ( ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) * ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_p] * F_inh * yout[ID_C_inh] ) + ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) * yout[ID_Q_cardiac] * yout[ID_Free] * R_ven_in - ( yout[ID_Q_p] * F_inh * yout[ID_C_inh] ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * pow ( yout[ID_Q_cardiac] * yout[ID_Free] , 2 ) - ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) ) ;

  yout[ID_C_pulv] = ( GE_ss ? C_pulv_ss : 0.0 ) ;

  yout[ID_C_art_in] = ( exist_lung ? yout[ID_C_lu] / P_lu : yout[ID_C_pulv] ) ;

  C_art_ss_lung = ( GE_ss ? ( ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_C_lu] ) / P_lu + yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * R_ven_in ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) : ( ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_C_lu] ) / P_lu + yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * R_ven_in ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) ) ;

  C_art_ss = ( exist_lung ? C_art_ss_lung : ( yout[ID_Q_cardiac] * yout[ID_Free] * ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_p] * F_inh * yout[ID_C_inh] ) + R_ven_in * pow ( yout[ID_Q_cardiac] * yout[ID_Free] , 2 ) + yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * R_ven_in ) / ( ( yout[ID_Q_cardiac] + Ven_loss_c ) * ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * ( yout[ID_Q_cardiac] + yout[ID_Free] * yout[ID_k_GFR] ) - yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * pow ( yout[ID_Q_cardiac] * yout[ID_Free] , 2 ) - ( yout[ID_Q_cardiac] * yout[ID_Free] + yout[ID_Q_p] / P_B ) * pow ( yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) , 2 ) ) ) ;

  yout[ID_C_art_comp] = ( arterial_ss ? C_art_ss : y[ID_A_art] / yout[ID_V_art] ) ;

  yout[ID_C_art] = ( single_blood ? y[ID_A_bl] / yout[ID_V_bl] : yout[ID_C_art_comp] ) ;

  yout[ID_R_glgi] = k_absgi * y[ID_A_glumen] ;

  yout[ID_R_glli] = k_absli * y[ID_A_glumen] / yout[ID_CarD] ;

  yout[ID_R_gili] = k_absli2 * y[ID_A_gi] / yout[ID_CarD] ;

  yout[ID_R_glfst] = k_unabs * y[ID_A_glumen] ;

  yout[ID_R_gifst] = k_fst * y[ID_A_gi] ;

  yout[ID_R_biliary] = yout[ID_Free] * yout[ID_k_bile] * y[ID_A_li] / P_li ;

  yout[ID_R_feces] = yout[ID_k_feces] * y[ID_A_fst] ;

  yout[ID_R_filust] = yout[ID_k_ust] * yout[ID_V_fil] * yout[ID_C_fil] ;

  R_venust = k_ven_ust * ( single_blood ? yout[ID_V_bl] * yout[ID_C_bl] : yout[ID_V_ven] * yout[ID_C_ven] ) ;

  yout[ID_R_urine] = yout[ID_k_urine] * y[ID_A_ust] ;

  yout[ID_R_artfil] = yout[ID_Free] * yout[ID_k_GFR] * yout[ID_C_art] ;

  yout[ID_R_reabs] = yout[ID_V_max_reabs] * yout[ID_C_fil] / ( K_m_reabs + yout[ID_C_fil] ) ;

  yout[ID_R_libind] = ( V_max_bind_li * yout[ID_C_li] * yout[ID_Free] / P_li ) / ( K_m_bind_li + yout[ID_C_li] * yout[ID_Free] / P_li ) ;

  yout[ID_R_liunbind] = k_off_li * y[ID_A_lib] ;

  yout[ID_R_blloss] = k_loss * ( single_blood ? yout[ID_C_bl] : yout[ID_C_ven] ) ;

  yout[ID_R_met_sat_lu] = y[ID_V_max_lu_t] * ( yout[ID_C_lu] / P_lu ) / ( K_m_lu + ( yout[ID_C_lu] / P_lu ) ) ;

  yout[ID_R_met_1st_lu] = yout[ID_k_met_lu] * ( yout[ID_C_lu] / P_lu ) * yout[ID_V_lu] ;

  yout[ID_R_met_sat_li] = y[ID_V_max_li_t] * ( yout[ID_C_li] / P_li ) / ( K_m_li + ( yout[ID_C_li] / P_li ) ) ;

  yout[ID_R_met_1st_li] = yout[ID_k_met_li] * ( yout[ID_C_li] / P_li ) * yout[ID_V_li] ;

  yout[ID_R_met_sat_om] = y[ID_V_max_om_t] * ( yout[ID_C_om] / P_om ) / ( K_m_om + ( yout[ID_C_om] / P_om ) ) ;

  yout[ID_R_met_1st_om] = yout[ID_k_met_om] * ( yout[ID_C_om] / P_om ) * yout[ID_V_om] ;

  R_enz_loss_lu = k_enz_loss * ( y[ID_V_max_lu_t] / yout[ID_V_lu] ) * yout[ID_R_met_sat_lu] ;

  R_enz_loss_li = k_enz_loss * ( y[ID_V_max_li_t] / yout[ID_V_li] ) * yout[ID_R_met_sat_li] ;

  R_enz_loss_om = k_enz_loss * ( y[ID_V_max_om_t] / yout[ID_V_om] ) * yout[ID_R_met_sat_om] ;

  yout[ID_R_inh] = yout[ID_Q_p] * F_inh * yout[ID_C_inh] ;

  yout[ID_R_exh] = yout[ID_Q_p] * ( GE_ss ? yout[ID_C_pulv] : yout[ID_C_lu] / P_lu ) / P_B ;

  R_ch_loss = KL * y[ID_A_ch] ;

  R_Qgi_in = yout[ID_Free] * yout[ID_Q_gi] * yout[ID_C_art] ;

  R_Qgi_out = yout[ID_Free] * yout[ID_Q_gi] * yout[ID_C_gi] / P_gi ;

  R_Qli_in = yout[ID_Free] * yout[ID_Q_li] * yout[ID_C_art] ;

  R_Qli_out = yout[ID_Free] * ( yout[ID_Q_li] + yout[ID_Q_gi] ) * yout[ID_C_li] / P_li ;

  R_Qki_in = yout[ID_Free] * yout[ID_Q_ki] * yout[ID_C_art] ;

  R_Qki_out = yout[ID_Free] * yout[ID_Q_ki] * yout[ID_C_ki] / P_ki ;

  R_Qtc1_in = yout[ID_Free] * yout[ID_Q_tc1] * yout[ID_C_art] ;

  R_Qtc1_out = yout[ID_Free] * yout[ID_Q_tc1] * yout[ID_C_tc1] / P_tc1 ;

  R_Qtc2_in = yout[ID_Free] * yout[ID_Q_tc2] * yout[ID_C_art] ;

  R_Qtc2_out = yout[ID_Free] * yout[ID_Q_tc2] * yout[ID_C_tc2] / P_tc2 ;

  R_Qtc3_in = yout[ID_Free] * yout[ID_Q_tc3] * yout[ID_C_art] ;

  R_Qtc3_out = yout[ID_Free] * yout[ID_Q_tc3] * yout[ID_C_tc3] / P_tc3 ;

  R_Qtc4_in = yout[ID_Free] * yout[ID_Q_tc4] * yout[ID_C_art] ;

  R_Qtc4_out = yout[ID_Free] * yout[ID_Q_tc4] * yout[ID_C_tc4] / P_tc4 ;

  R_Qtc5_in = yout[ID_Free] * yout[ID_Q_tc5] * yout[ID_C_art] ;

  R_Qtc5_out = yout[ID_Free] * yout[ID_Q_tc5] * yout[ID_C_tc5] / P_tc5 ;

  R_Qom_in = yout[ID_Free] * yout[ID_Q_om] * yout[ID_C_art] ;

  R_Qom_out = yout[ID_Free] * yout[ID_Q_om] * yout[ID_C_om] / P_om ;

  R_Qrb_in = yout[ID_Free] * yout[ID_Q_rb] * yout[ID_C_art] ;

  R_Qrb_out = yout[ID_Free] * yout[ID_Q_rb] * yout[ID_C_rb] / P_rb ;

  yout[ID_A_body] = y[ID_A_glumen] + y[ID_A_bl] + y[ID_A_ven] + y[ID_A_art] + y[ID_A_lu] + y[ID_A_gi] + y[ID_A_li] + y[ID_A_ki] + y[ID_A_fil] + y[ID_A_ust] + y[ID_A_fst] + y[ID_A_tc1] + y[ID_A_tc2] + y[ID_A_tc3] + y[ID_A_tc4] + y[ID_A_tc5] + y[ID_A_om] + y[ID_A_rb] + y[ID_A_lib] ;

  yout[ID_A_out] = y[ID_A_urine] + y[ID_A_fecal] + y[ID_A_loss] + y[ID_A_exh] + y[ID_A_met_sat] + y[ID_A_met_1st] ;

  yout[ID_A_bal] = y[ID_A_in] - yout[ID_A_body] - yout[ID_A_out] ;

  Q_body = yout[ID_Q_gi] + yout[ID_Q_li] + yout[ID_Q_ki] + yout[ID_Q_tc1] + yout[ID_Q_tc2] + yout[ID_Q_tc3] + yout[ID_Q_tc4] + yout[ID_Q_tc5] + yout[ID_Q_om] + yout[ID_Q_rb] ;

  yout[ID_Q_bal] = yout[ID_Q_cardiac] - Q_body ;

  yout[ID_Q_relative] = yout[ID_Q_cardiac] / Q_body - 1 ;

  yout[ID_V_total_bl_value] = ( single_blood ? yout[ID_V_bl] : ( ( 1 - venous_ss ) * yout[ID_V_ven] + ( 1 - arterial_ss ) * yout[ID_V_art] ) ) ;

  yout[ID_V_lu_value] = ( exist_lung ? yout[ID_V_lu] : 0 ) ;

  yout[ID_V_ki_value] = ( Q_kic ? yout[ID_V_ki] : 0 ) ;

  yout[ID_V_fil_value] = ( yout[ID_k_GFR] ? yout[ID_V_fil] : 0 ) ;

  yout[ID_V_li_value] = ( Q_lic ? yout[ID_V_li] : 0 ) ;

  yout[ID_V_gi_value] = ( Q_gic ? yout[ID_V_gi] : 0 ) ;

  yout[ID_V_tc1_value] = ( Q_tc1c ? yout[ID_V_tc1] : 0 ) ;

  yout[ID_V_tc2_value] = ( Q_tc2c ? yout[ID_V_tc2] : 0 ) ;

  yout[ID_V_tc3_value] = ( Q_tc3c ? yout[ID_V_tc3] : 0 ) ;

  yout[ID_V_tc4_value] = ( Q_tc4c ? yout[ID_V_tc4] : 0 ) ;

  yout[ID_V_tc5_value] = ( Q_tc5c ? yout[ID_V_tc5] : 0 ) ;

  yout[ID_V_om_value] = ( Q_omc ? yout[ID_V_om] : 0 ) ;

  yout[ID_V_rb_value] = ( Q_rbc ? yout[ID_V_rb] : 0 ) ;

  yout[ID_total_volume_frac] = ( yout[ID_V_total_bl_value] + yout[ID_V_lu_value] + yout[ID_V_ki_value] + yout[ID_V_fil_value] + yout[ID_V_li_value] + yout[ID_V_gi_value] + yout[ID_V_tc1_value] + yout[ID_V_tc2_value] + yout[ID_V_tc3_value] + yout[ID_V_tc4_value] + yout[ID_V_tc5_value] + yout[ID_V_om_value] + yout[ID_V_rb_value] ) / yout[ID_BW_out] ;

  ydot[ID_A_bl] = ( single_blood ? R_Qli_out + R_Qki_out + R_Qtc1_out + R_Qtc2_out + R_Qtc3_out + R_Qtc4_out + R_Qtc5_out + R_Qom_out + R_Qrb_out - R_Qli_in - R_Qki_in - R_Qgi_in - R_Qtc1_in - R_Qtc2_in - R_Qtc3_in - R_Qtc4_in - R_Qtc5_in - R_Qom_in - R_Qrb_in - yout[ID_R_artfil] + y[ID_R_IV] - R_venust - yout[ID_R_blloss] : 0.0 ) ;

  ydot[ID_AUC_bl] = yout[ID_C_bl] ;

  ydot[ID_A_loss] = yout[ID_R_blloss] ;

  ydot[ID_A_ven] = ( 1 - single_blood ) * ( venous_ss ? 0.0 : R_Qli_out + R_Qki_out + R_Qtc1_out + R_Qtc2_out + R_Qtc3_out + R_Qtc4_out + R_Qtc5_out + R_Qom_out + R_Qrb_out + y[ID_R_IV] + yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * yout[ID_C_art] - yout[ID_Q_cardiac] * yout[ID_C_ven] - R_venust - yout[ID_R_blloss] ) ;

  ydot[ID_AUC_ven] = yout[ID_C_ven] ;

  ydot[ID_A_art] = ( 1 - single_blood ) * ( arterial_ss ? 0.0 : yout[ID_Q_cardiac] * yout[ID_Free] * yout[ID_C_art_in] + yout[ID_Q_cardiac] * ( 1 - yout[ID_Free] ) * yout[ID_C_ven] - yout[ID_Q_cardiac] * yout[ID_C_art] - yout[ID_R_artfil] ) ;

  ydot[ID_AUC_art] = yout[ID_C_art] ;

  ydot[ID_Conc] = 0 ;

  ydot[ID_Q_cc] = 0 ;

  ydot[ID_A_ch] = ( V_ch ? NCH * ( yout[ID_R_exh] - yout[ID_R_inh] ) - R_ch_loss : 0.0 ) ;

  ydot[ID_A_exh] = yout[ID_R_exh] ;

  ydot[ID_A_lu] = exist_lung * ( ( GE_ss ? yout[ID_Q_cardiac] * yout[ID_Free] * ( yout[ID_C_pulv] - yout[ID_C_lu] / P_lu ) : yout[ID_R_inh] - yout[ID_R_exh] + yout[ID_Q_cardiac] * yout[ID_Free] * ( yout[ID_C_ven] - yout[ID_C_lu] / P_lu ) ) - yout[ID_R_met_sat_lu] - yout[ID_R_met_1st_lu] ) ;

  ydot[ID_AUC_lu] = yout[ID_C_lu] ;

  ydot[ID_A_glumen] = y[ID_R_oral] - yout[ID_R_glgi] - yout[ID_R_glli] - yout[ID_R_glfst] ;

  ydot[ID_A_gi] = R_Qgi_in - R_Qgi_out - yout[ID_R_gifst] + yout[ID_R_glgi] - yout[ID_R_gili] ;

  ydot[ID_AUC_gi] = yout[ID_C_gi] ;

  ydot[ID_A_li] = R_Qli_in + R_Qgi_out - R_Qli_out - yout[ID_R_biliary] + yout[ID_R_glli] + yout[ID_R_gili] + R_0bgli - yout[ID_R_met_sat_li] - yout[ID_R_met_1st_li] - yout[ID_R_libind] + yout[ID_R_liunbind] ;

  ydot[ID_AUC_li] = yout[ID_C_li] ;

  ydot[ID_A_fst] = yout[ID_R_biliary] + yout[ID_R_glfst] + yout[ID_R_gifst] - yout[ID_R_feces] ;

  ydot[ID_A_fecal] = yout[ID_R_feces] ;

  ydot[ID_A_lib] = yout[ID_R_libind] - yout[ID_R_liunbind] ;

  ydot[ID_A_ki] = R_Qki_in - R_Qki_out + yout[ID_R_reabs] ;

  ydot[ID_AUC_ki] = yout[ID_C_ki] ;

  ydot[ID_A_fil] = yout[ID_R_artfil] - yout[ID_R_reabs] - yout[ID_R_filust] ;

  ydot[ID_A_ust] = yout[ID_R_filust] + R_venust - yout[ID_R_urine] ;

  ydot[ID_A_urine] = yout[ID_R_urine] ;

  ydot[ID_A_tc1] = R_Qtc1_in - R_Qtc1_out ;

  ydot[ID_AUC_tc1] = yout[ID_C_tc1] ;

  ydot[ID_A_tc2] = R_Qtc2_in - R_Qtc2_out ;

  ydot[ID_AUC_tc2] = yout[ID_C_tc2] ;

  ydot[ID_A_tc3] = R_Qtc3_in - R_Qtc3_out ;

  ydot[ID_AUC_tc3] = yout[ID_C_tc3] ;

  ydot[ID_A_tc4] = R_Qtc4_in - R_Qtc4_out ;

  ydot[ID_AUC_tc4] = yout[ID_C_tc4] ;

  ydot[ID_A_tc5] = R_Qtc5_in - R_Qtc5_out ;

  ydot[ID_AUC_tc5] = yout[ID_C_tc5] ;

  ydot[ID_A_om] = R_Qom_in - R_Qom_out - yout[ID_R_met_sat_om] - yout[ID_R_met_1st_om] ;

  ydot[ID_AUC_om] = yout[ID_C_om] ;

  ydot[ID_V_max_lu_t] = k_enz_resyn * ( yout[ID_V_max_lu] - y[ID_V_max_lu_t] ) - R_enz_loss_lu ;

  ydot[ID_V_max_li_t] = k_enz_resyn * ( yout[ID_V_max_li] - y[ID_V_max_li_t] ) - R_enz_loss_li ;

  ydot[ID_V_max_om_t] = k_enz_resyn * ( yout[ID_V_max_om] - y[ID_V_max_om_t] ) - R_enz_loss_om ;

  ydot[ID_A_rb] = R_Qrb_in - R_Qrb_out ;

  ydot[ID_AUC_rb] = yout[ID_C_rb] ;

  ydot[ID_A_met_sat_li] = yout[ID_R_met_sat_li] ;

  ydot[ID_A_met_sat_om] = yout[ID_R_met_sat_om] ;

  ydot[ID_A_met_sat_lu] = yout[ID_R_met_sat_lu] ;

  ydot[ID_A_met_sat] = yout[ID_R_met_sat_li] + yout[ID_R_met_sat_om] + yout[ID_R_met_sat_lu] ;

  ydot[ID_A_met_1st_li] = yout[ID_R_met_1st_li] ;

  ydot[ID_A_met_1st_lu] = yout[ID_R_met_1st_lu] ;

  ydot[ID_A_met_1st_om] = yout[ID_R_met_1st_om] ;

  ydot[ID_A_met_1st] = yout[ID_R_met_1st_li] + yout[ID_R_met_1st_lu] + yout[ID_R_met_1st_om] ;

  ydot[ID_R_IV] = 0 ;

  ydot[ID_R_oral] = 0 ;

  ydot[ID_A_in] = y[ID_R_IV] + y[ID_R_oral] + yout[ID_R_inh] + R_0bgli ;

} /* derivs */


/*----- Jacobian calculations: */
void jac (int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{

} /* jac */


/*----- Events calculations: */
void event (int *n, double *t, double *y)
{

} /* event */

/*----- Roots calculations: */
void root (int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip)
{

} /* root */

