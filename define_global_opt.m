function [] = define_global_opt()
% This function is used to define global variables used in the simulation.

global iPO2 idPsi iH_x iATP_x iADP_x iAMP_x iGTP_x iGDP_x iPI_x iNADH_x iQH2_x ...
       iOAA_x iACCOA_x iCIT_x iICIT_x iAKG_x iSCOA_x iCOASH_x iSUC_x iFUM_x ...
       iMAL_x iGLU_x iASP_x iK_x iMg_x iCO2tot_x iCred_i iATP_i iADP_i iAMP_i ...
       iPI_i iH_i iMg_i iK_i iATP_c iADP_c iPI_c iH_c iMg_c iK_c iPYR_x iPYR_i ...
       iPYR_c iCIT_i iCIT_c iAKG_i iAKG_c iSUC_i iSUC_c iMAL_i iMAL_c iASP_i ...
       iASP_c iGLU_i iGLU_c iFUM_i iFUM_c iICIT_i iICIT_c iPCr_c iAMP_c iCr_c;
     
global i_J_C1 i_J_C3 i_J_C4 i_J_F1 i_J_ANT i_J_PI1 i_J_Hle i_J_KH i_J_pdh i_J_cits ...
       i_J_acon i_J_isod i_J_akgd i_J_scoas i_J_sdh i_J_fum i_J_mdh i_J_ndk i_J_got ...
	     i_J_PI2 i_J_ATPt i_J_ADPt i_J_AMPt i_J_PYR_H i_J_GLU_H i_J_CIT_MAL i_J_SUC_MAL ... 
       i_J_AKG_MAL i_J_MAL_PI i_J_ASP_GLU i_J_PYRt i_J_CITt i_J_ICITt i_J_MALt ... 
	     i_J_AKGt i_J_SUCt i_J_FUMt i_J_GLUt i_J_ASPt i_J_CKe i_J_AKe i_J_AtC;
     
global iC16Carn_m iC16CoA_m iC16EnoylCoA_m iC16OHCoA_m iC16KetoCoA_m...   
       iC14Carn_cy iC14Carn_m iC14CoA_m iC14EnoylCoA_m iC14OHCoA_m...    
       iC14KetoCoA_m iC12Carn_cy iC12Carn_m iC12CoA_m iC12EnoylCoA_m iC12OHCoA_m...   
       iC12KetoCoA_m iC10Carn_cy iC10Carn_m iC10CoA_m iC10EnoylCoA_m iC10OHCoA_m...  
       iC10KetoCoA_m iC8Carn_cy iC8Carn_m iC8CoA_m iC8EnoylCoA_m iC8OHCoA_m...  
       iC8KetoCoA_m iC6Carn_cy iC6Carn_m iC6CoA_m iC6EnoylCoA_m iC6OHCoA_m...
       iC6KetoCoA_m iC4Carn_cy iC4Carn_m iC4CoA_m...
       iC4EnoylCoA_m iC4OHCoA_m iC4KetoCoA_m iAcetylCoAMAT iFADH_m iNADHm iCoAMAT...
       iC16AcylCoACYT iC16Carn_cy;

%% define indices for all state variables
%   For biochemical concentration variables the notations
%   "_x", "_i", and "_c" denote matrix, im space, and cytoplasm,
%   respectively. Thus iATP_x indexes ATP in the matrix.
%   
% (i) oxygen species
iPO2 = 1;
% (ii) Matrix species and dPsi
idPsi        = 2;
iH_x         = 3; 
iATP_x       = 4;
iADP_x       = 5;
iAMP_x       = 6;
iGTP_x       = 7;
iGDP_x       = 8;
iPI_x        = 9;
iNADH_x      = 10;
iQH2_x       = 11;
iOAA_x       = 12;
iACCOA_x     = 13;
iCIT_x       = 14;
iICIT_x      = 15;
iAKG_x       = 16;
iSCOA_x      = 17;
iCOASH_x     = 18;
iSUC_x       = 19;
iFUM_x       = 20;
iMAL_x       = 21;
iGLU_x       = 22;
iASP_x       = 23;
iK_x         = 24;
iMg_x        = 25;
iCO2tot_x    = 26;
%  (iii) IM space species
iCred_i      = 27;
iATP_i       = 28;
iADP_i       = 29;
iAMP_i       = 30;
iPI_i        = 31;
iH_i         = 32;
iMg_i        = 33;
iK_i         = 34;
%  (iv) Cytoplasmic species
iATP_c       = 35;
iADP_c       = 36;
iPI_c        = 37;
iH_c         = 38;
iMg_c        = 39;
iK_c         = 40;
% others
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
iPCr_c      = 60;
iAMP_c      = 61;
iCr_c = 62;

iC16Carn_cy     = 63;    % C16 AcylCarnitine cytosol
iC16Carn_m      = 64;    % C16 AcylCarnitine Matrix
iC16CoA_m       = 65;    % C16 AcylCoA Matrix
iC16EnoylCoA_m  = 66;    % C16 EnoylCoA Matrix
iC16OHCoA_m     = 67;    % C16 HydroxoxyacylCoA Matrix
iC16KetoCoA_m   = 68;    % C16 KetoacylCoA Matrix

iC14Carn_cy     = 69;    % C14 AcylCarnitine cytosol
iC14Carn_m      = 70;    % C14 AcylCarnitine Matrix
iC14CoA_m       = 71;    % C14 AcylCoA Matrix
iC14EnoylCoA_m  = 72;    % C14 EnoylCoA Matrix
iC14OHCoA_m     = 73;    % C14 HydroxoxyacylCoA Matrix
iC14KetoCoA_m   = 74;

iC12Carn_cy     = 75;  
iC12Carn_m      = 76;    
iC12CoA_m       = 77;
iC12EnoylCoA_m  = 78;    
iC12OHCoA_m     = 79;    
iC12KetoCoA_m   = 80; 

iC10Carn_cy     = 81; 
iC10Carn_m      = 82;   
iC10CoA_m       = 83;
iC10EnoylCoA_m  = 84;
iC10OHCoA_m     = 85;
iC10KetoCoA_m   = 86; 

iC8Carn_cy      = 87; 
iC8Carn_m       = 88;    
iC8CoA_m        = 89;
iC8EnoylCoA_m   = 90;
iC8OHCoA_m      = 91;
iC8KetoCoA_m    = 92; 

iC6Carn_cy      = 93; 
iC6Carn_m       = 94;    
iC6CoA_m        = 95;
iC6EnoylCoA_m   = 96;
iC6OHCoA_m      = 97;
iC6KetoCoA_m    = 98; 

iC4Carn_cy      = 99; 
iC4Carn_m       = 100;     
iC4CoA_m        = 101;
iC4EnoylCoA_m   = 102;
iC4OHCoA_m      = 103;
iC4KetoCoA_m    = 104; 

%State variables that are involved in both TCA/Oxphos and b oxidation
iAcetylCoAMAT   = 105; %
iFADH_m         = 106;
iNADHm          = 107; 
iCoAMAT         = 108; 
iC16AcylCoACYT  = 109;


%%  Defining indices for all fluxes
	i_J_C1 = 4 ;
	i_J_C3 = 5;
	i_J_C4 = 6;
	i_J_F1 = 7; 
	i_J_ANT = 8; 
	i_J_PI1 = 9;
	i_J_Hle = 10; 
	i_J_KH = 11;
	i_J_pdh = 12;
	i_J_cits = 13;
    i_J_acon = 14;
	i_J_isod = 15; 
	i_J_akgd = 16;
	i_J_scoas = 17;
	i_J_sdh = 18;
	i_J_fum = 19; 
	i_J_mdh = 20;
	i_J_ndk = 21;
	i_J_got = 22;
	i_J_PI2 = 23; 
	i_J_ATPt = 24;
	i_J_ADPt = 25; 
	i_J_AMPt = 26; 
    i_J_PYR_H = 27; 
	i_J_GLU_H = 28;
	i_J_CIT_MAL = 29; 
	i_J_SUC_MAL = 30; 
	i_J_AKG_MAL = 31; 
	i_J_MAL_PI = 32;	
	i_J_ASP_GLU = 33; 
	i_J_PYRt = 34; 
	i_J_CITt = 35; 
	i_J_ICITt = 36; 
	i_J_MALt = 37; 
	i_J_AKGt = 38; 
	i_J_SUCt = 39;
	i_J_FUMt = 40; 
	i_J_GLUt = 41;
	i_J_ASPt = 42;
	i_J_CKe = 43;
	i_J_AKe = 44;
	i_J_AtC = 45;