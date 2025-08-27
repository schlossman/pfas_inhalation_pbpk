initParms <- function(newParms = NULL) {
  parms <- c(
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
    V_max_lu_0 = 0.0
  )

  if (!is.null(newParms)) {
    if (!all(names(newParms) %in% c(names(parms)))) {
      stop("illegal parameter name")
    }
    parms[names(newParms)] <- newParms
  }

  parms <- within(as.list(parms), {
  })
  out <- .C("getParms",  as.double(parms),
            out=double(length(parms)),
            as.integer(length(parms)))$out
  names(out) <- names(parms)
  out
}

Outputs <- c(
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
    "total_volume_frac"
)

initStates <- function(parms, newStates = NULL) {
  Y <- c(
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
    Q_cc = 0.0
  )

  Y <- within(c(as.list(parms),as.list(Y)), {
    Y["A_bl"] <- A_bl0 
    Y["A_ven"] <- A_ven0 
    Y["A_art"] <- 0.0 
    Y["A_lu"] <- 0.0 
    Y["A_gi"] <- 0.0 
    Y["A_li"] <- 0.0 
    Y["A_fst"] <- oral_dose_init * BW * F_unabs 
    Y["A_ki"] <- 0.0 
    Y["A_fil"] <- 0.0 
    Y["A_ust"] <- 0.0 
    Y["A_tc1"] <- 0.0 
    Y["A_tc2"] <- 0.0 
    Y["A_tc3"] <- 0.0 
    Y["A_tc4"] <- 0.0 
    Y["A_tc5"] <- 0.0 
    Y["A_om"] <- 0.0 
    Y["A_rb"] <- 0.0 
    Y["A_glumen"] <- oral_dose_init * BW * ( 1 - F_unabs ) 
    Y["R_IV"] <- 0.0 
    Y["R_oral"] <- 0.0 
    Y["A_loss"] <- 0.0 
    Y["A_exh"] <- 0.0 
    Y["A_urine"] <- 0.0 
    Y["A_fecal"] <- 0.0 
    Y["A_in"] <- oral_dose_init * BW + iv_dose * BW 
    Y["A_lib"] <- 0.0 
    Y["Conc"] <- Conc_ambient 
    Y["Q_cc"] <- Q_cardiacc 
    Y["A_ch"] <- Conc_init * GASD * V_ch 

    Y["AUC_bl"] <- 0.0 
    Y["AUC_ven"] <- 0.0 
    Y["AUC_art"] <- 0.0 
    Y["AUC_lu"] <- 0.0 
    Y["AUC_gi"] <- 0.0 
    Y["AUC_li"] <- 0.0 
    Y["AUC_ki"] <- 0.0 
    Y["AUC_tc1"] <- 0.0 
    Y["AUC_tc2"] <- 0.0 
    Y["AUC_tc3"] <- 0.0 
    Y["AUC_tc4"] <- 0.0 
    Y["AUC_tc5"] <- 0.0 
    Y["AUC_om"] <- 0.0 
    Y["AUC_rb"] <- 0.0 

    Y["A_met_sat"] <- 0.0 
    Y["A_met_sat_li"] <- 0.0 
    Y["A_met_sat_om"] <- 0.0 
    Y["A_met_sat_lu"] <- 0.0 
    Y["A_met_1st"] <- 0.0 
    Y["A_met_1st_li"] <- 0.0 
    Y["A_met_1st_lu"] <- 0.0 
    Y["A_met_1st_om"] <- 0.0 

    Y["V_max_om_t"] <- V_max_om_0 
    Y["V_max_li_t"] <- V_max_li_0 
    Y["V_max_lu_t"] <- V_max_lu_0 

  })$Y

  if (!is.null(newStates)) {
    if (!all(names(newStates) %in% c(names(Y)))) {
      stop("illegal state variable name in newStates")
    }
    Y[names(newStates)] <- newStates
  }

.C("initState", as.double(Y));
Y
}
