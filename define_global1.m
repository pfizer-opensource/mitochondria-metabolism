function define_global()

% Define indices of model variables as global variables

global iH_x idPsi iATP_x iADP_x iAMP_x iGTP_x iGDP_x iPi_x iNADH_x ...
       iQH2_x iOAA_x iACCOA_x iCIT_x iICIT_x iAKG_x iSCOA_x iCOASH_x ...
       iSUC_x iFUM_x iMAL_x iGLU_x iASP_x iK_x iMg_x iO2_x iCO2tot_x ...
       iCred_i iATP_i iADP_i iAMP_i iPi_i iH_i iMg_i iK_i ...
       iATP_c iADP_c iPi_c iH_c iMg_c iK_c ...
       iPYR_x iPYR_i iPYR_c iCIT_i iCIT_c iAKG_i iAKG_c iSUC_i iSUC_c ...
       iMAL_i iMAL_c iASP_i iASP_c iGLU_i iGLU_c iFUM_i iFUM_c ...
       iICIT_i iICIT_c iGLC_c iG6P_c iPCr_c iAMP_c;


%% Define indices for all state variables
%  (i) Matrix species and dPsi
iH_x         = 1; 
idPsi        = 2;
iATP_x       = 3;
iADP_x       = 4;
iAMP_x       = 5;
iGTP_x       = 6;
iGDP_x       = 7;
iPi_x        = 8;
iNADH_x      = 9;
iQH2_x       = 10;
iOAA_x       = 11;
iACCOA_x     = 12;
iCIT_x       = 13;
iICIT_x      = 14;
iAKG_x       = 15;
iSCOA_x      = 16;
iCOASH_x     = 17;
iSUC_x       = 18;
iFUM_x       = 19;
iMAL_x       = 20;
iGLU_x       = 21;
iASP_x       = 22;
iK_x         = 23;
iMg_x        = 24;
iO2_x        = 25;
iCO2tot_x    = 26;
%  (ii) IM space species
iCred_i      = 27;
iATP_i       = 28;
iADP_i       = 29;
iAMP_i       = 30;
iPi_i        = 31;
iH_i         = 32;
iMg_i        = 33;
iK_i         = 34;
%  (iii) Cytoplasmic species
iATP_c       = 35;
iADP_c       = 36;
iPi_c        = 37;
iH_c         = 38;
iMg_c        = 39;
iK_c         = 40;

% added variable
iPYR_x      = 41;
iPYR_i      = 42;
iPYR_c      = 43;
iCIT_i      = 44;
iCIT_c      = 45;
iAKG_i      = 46;
iAKG_c      = 47;
iSUC_i      = 48;
iSUC_c      = 49;
iMAL_i      = 50;
iMAL_c      = 51;
iASP_i      = 52;
iASP_c      = 53;
iGLU_i      = 54;
iGLU_c      = 55;
iFUM_i      = 56;
iFUM_c      = 57;
iICIT_i     = 58;
iICIT_c     = 59;
iGLC_c      = 60;
iG6P_c      = 61;
iPCr_c      = 62;
iAMP_c      = 63;
