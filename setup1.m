% This setup program is used to load experimental data, certain model 
% parameters, and initial conditions for simulation.


%clear;
format short e;
define_global_opt; % This script defines global variables.

% %% Loading experimental data
% % These are experimental data from physiological hearts under different
% % work rates. These data are collected from a series of publications from
% % Zhang Lab at University of Minnesota. 
% % Please refer to Wu et al., J Physiol 586 (17), 4193 (2008) for detailed 
% % data sources.
% load HEP_MVO2.txt;
% MVO2_st_z = HEP_MVO2(:,1); % average values of steady-state oxygen consumption rates
% MVO2_st_e = HEP_MVO2(:,2); % standard deviations of steady-state oxygen consumption rates
% PCrATP_st_z = HEP_MVO2(:,3); % average values of steady-state phosphocreatine-to-ATP ratios
% PCrATP_st_e = HEP_MVO2(:,4); % standard deviations of steady-state phosphocreatine-to-ATP ratios
% PIPCr_st_z = HEP_MVO2(:,5); % average values of steady-state Pi-to-PCr ratios
% PIPCr_st_e = HEP_MVO2(:,6); % standard deviations of steady-state Pi-to-PCr ratios 
% 
% % These are experimental data from left ventricular hypertropy (LVH) hearts 
% % under different work rates. These data are collected from a series of 
% % publications from Bache Lab at University of Minnesota. 
% % Please refer to Wu et al., PNAS 106 (17), 7143 (2009) for detailed data 
% % sources. 
% load HEP_MVO2_LVH.txt;
% % HEP_MVO2_LVH (1:4,:) are the data on normal hearts, Bache et al., Cardio Res, 1999
% % HEP_MVO2_LVH (5:8,:) are the data on LVH hearts, Bache et al, Cardio Res, 1999
% % HEP_MVO2_LVH (9:11,:) are the data on LVH hearts, Bache et al, AJP-Heart, 1994
% MVO2_st_LVH_z = HEP_MVO2_LVH(5:11,1);
% MVO2_st_LVH_e = HEP_MVO2_LVH(5:11,2);
% PCrATP_st_LVH_z = HEP_MVO2_LVH(5:11,3);
% PCrATP_st_LVH_e = HEP_MVO2_LVH(5:11,4);
% PIPCr_st_LVH_z = HEP_MVO2_LVH(5:11,5);
% PIPCr_st_LVH_e = HEP_MVO2_LVH(5:11,6);


%% Setting up myocardial density and composition 
% Obtained from Vinnakota and Bassingthwaighte,AJP-Heart Circ Physiol 286 
% (5) H1742 (2004).

% Density, g region/ml region
rou_water = 1.0; % assumed water density
rou_tissue = 1.053; 
rou_plasma = 1.023; 
rou_RBC = 1.097; 
rou_ISF = 1.022; 
rou_cell = 1.06; 
rou_cyto = 1.044; 
rou_mito = 1.09;
rou_SR = 1.137; 
% Volume, ml region/g tissue
Vtissue_tissue = 0.95;
Vplasma_tissue = 0.056;
VRBC_tissue = 0.033;
VISF_tissue = 0.167;
Vcell_tissue = 0.694;
Vcyto_tissue = 0.472;
Vmito_tissue = 0.2;
VSR_tissue = 0.024;
% Water fraction, g water/g region
W_tissue = 0.79; 
W_plasma = 0.924;
W_RBC = 0.67;
W_ISF = 0.924;
W_cell = 0.755;
W_cyto = 0.807;
W_mito = 0.664;
W_SR = 0.518;
% Protein fraction, g protein/g region
P_tissue = 0.176;
P_plasma = 0.059;
P_RBC = 0.308;
P_ISF = 0.065;
P_cell = 0.204;
P_cyto = 0.17;
P_mito = 0.25;
P_SR = 0.445;

% ---------------------------------------------------
Vmito = 0.2882;                 % (ml mito / ml cell)
Vcyto = 0.6801;                 % (ml cyto / ml cell)
Rm_cyto = Vmito / Vcyto;        % Volume ratio mito volume / cytoplasm volume
Rm_cell = Vmito;                % Volume ratio mito volume / cell volume
Rc_cell = Vcyto;                % Volume ratio of cytosol / cell volume
W_c = 0.807*1.044;              % cytosol water space (ml water per ml cytosol) [VB2002]
W_m = 0.664*1.09;               % mitochondrial water space (ml water per ml mito) [VB2002]
W_x = 0.9*W_m;                  % Matrix water space
W_i = 0.1*W_m;                  % IM water space

Rc_t = 0.73078; % ratio volume cell/volume tissue
rho = 1.05; % tissue density grams per ml

% ---------------------------------------------------
%  (iii) Pooled concentrations
Ctot   = 2.70e-3;               % M; total cytoC
Qtot   = 1.35e-3;               % M; total Q + QH2
NADtot = 2.97e-3;               % M; total NAD+NADH
FADtot = 0.1e-3;              % M; total FADred+FADox (from Feng Yang)

%% Setting the initial concentrations to be the baseline ones obtained 
% Obtained from the control steady-state simulations at the baseline work
% rate

xo_con = [ 2.5000e+001
  1.7787e+002
  5.7261e-008
  2.8369e-003
  7.1631e-003
  1.0000e-006
  1.4340e-003
  3.5661e-003
  3.9155e-004
  2.4876e-003
  5.1367e-004
  5.1194e-011
  2.9782e-003
  4.2138e-004
  1.1630e-005
  1.8011e-007
  2.1823e-007
  9.4686e-006
  6.3456e-005
  7.4093e-005
  2.3412e-004
  1.0221e-004
  1.5307e-008
  9.3713e-002
  9.8624e-004
  2.1400e-002
  3.4883e-004
  1.1284e-002
  7.2404e-005
  4.0267e-007
  2.2800e-004
  7.9433e-008
  1.0030e-003
  1.3000e-001
  1.1281e-002
  7.4857e-005
  2.2863e-004
  7.9433e-008
  1.0030e-003
  1.3000e-001
  9.6854e-005
  7.4771e-005
  7.5000e-005
  1.6586e-004
  1.6586e-004
  9.7514e-008
  9.7514e-008
  3.5750e-005
  3.5750e-005
  1.2676e-004
  1.2676e-004
  1.1928e-005
  1.1928e-005
  7.3773e-005
  7.3773e-005
  1.5000e-005
  1.5000e-005
  1.5000e-005
  1.5000e-005
  2.6760e-002
  4.0267e-007
  3.4392e-002];

% %% Translating biopsy data into local concentrations
% % --------------------------------------------
% % Normal hearts: 
% CT_Biopsy = 115.8; % umol (g dry wt)^{-1}, Zhang et al., JCI, 1993
% CT_Biopsy_e = 3.0; % umol (g dry wt)^{-1}, Zhang et al., JCI, 1993
% ATPBL_Biopsy = 23.9; % umol (g dry wt)^{-1}, Zhang et al., JCI, 1993
% ATPBL_Biopsy_e = 0.84; % umol (g dry wt)^{-1}, Zhang et al., JCI, 1993
% CrPATPBL_cellwater = 2.119; % (mol (cell water)^{-1})/(mol (cell water)^{-1}), Wu et al., J Physiol, 2008
% CrPATPBL_cell = 2.119; % (mol (cell water)^{-1})/(mol (cell water)^{-1}), Wu et al., J Physiol, 2008
% CrPATPBL_cell_e = 0.108; % computed from the original data;
% 
% x_AtC = 0.36e-3; % Basal ATP consumption rate, mol s^{-1} (l cell)^{-1}, Wu et al., J Physiol, 2008
% x_AtC_BL = 0.36e-3; % Basal ATP consumption rate, mol s^{-1} (l cell)^{-1}, Wu et al., J Physiol, 2008
% x_AtC_max = 1.20e-3; % Basal ATP consumption rate, mol s^{-1} (l cell)^{-1}, Wu et al., J Physiol, 2008
%  
% CT_cytowater = CT_Biopsy*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% CT_cytowater_e = CT_Biopsy_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater = ATPBL_Biopsy*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater_e = ATPBL_Biopsy_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% 
% CT_cellwater = CT_Biopsy*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% CT_cellwater_e = CT_Biopsy_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater = ATPBL_Biopsy*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater_e = ATPBL_Biopsy_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% 
% CT_cell = CT_Biopsy*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell)^{-1}
% CT_cell_e = CT_Biopsy_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell)^{-1}
% ATPBL_cell = ATPBL_Biopsy*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell)^{-1}
% ATPBL_cell_e = ATPBL_Biopsy_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell)^{-1}
% 
% CrPBL_cell = CrPATPBL_cell*ATPBL_cell; % mol (l cell water)^{-1}
% CrPBL_cellwater = CrPATPBL_cellwater*ATPBL_cellwater; % mol (l cell water)^{-1}
% CrPBL_cytowater = CrPBL_cellwater*W_cell*rou_cell*Vcell_tissue/W_cyto/rou_cyto/Vcyto_tissue; % mol (l cyto water)^{-1}
% 
% % --------------------------------------------
% % Moderate LVH: (LVH hearts, Bache et al., AJP, 1994)
% CT_Biopsy_mLVH = 117.5; % umol (g dry wt)^{-1}
% ATPBL_Biopsy_mLVH = 22.5; % umol (g dry wt)^{-1}
% CT_Biopsy_mLVH_e = 3.3; % umol (g dry wt)^{-1}
% ATPBL_Biopsy_mLVH_e = 2.2; % umol (g dry wt)^{-1}
% CrPATPBL_cellwater_mLVH = 2.21; % Bache et al., ATP, 1994
% CrPATPBL_cell_mLVH = 2.21; % Bache et al., ATP, 1994
% CrPATPBL_cell_mLVH_e = 0.12; % Bache et al., ATP, 1994
% x_AtC_mLVH = 0.485e-3; % Corrected basal ATP consumption rate, mol s^{-1} (l cell)^{-1}
% 
% CT_cytowater_mLVH = CT_Biopsy_mLVH*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater_mLVH = ATPBL_Biopsy_mLVH*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% 
% CT_cellwater_mLVH = CT_Biopsy_mLVH*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater_mLVH = ATPBL_Biopsy_mLVH*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% 
% CT_cell_mLVH = CT_Biopsy_mLVH*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% ATPBL_cell_mLVH = ATPBL_Biopsy_mLVH*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% 
% CT_cytowater_mLVH_e = CT_Biopsy_mLVH_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater_mLVH_e = ATPBL_Biopsy_mLVH_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% CT_cellwater_mLVH_e = CT_Biopsy_mLVH_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater_mLVH_e = ATPBL_Biopsy_mLVH_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% CT_cell_mLVH_e = CT_Biopsy_mLVH_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% ATPBL_cell_mLVH_e = ATPBL_Biopsy_mLVH_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% 
% CrPBL_cell_mLVH = CrPATPBL_cell_mLVH*ATPBL_cell_mLVH; % mol (l cell water)^{-1}
% CrPBL_cellwater_mLVH = CrPATPBL_cellwater_mLVH*ATPBL_cellwater_mLVH; % mol (l cell water)^{-1}
% CrPBL_cytowater_mLVH = CrPBL_cellwater_mLVH*W_cell*rou_cell*Vcell_tissue/W_cyto/rou_cyto/Vcyto_tissue; % mol (l cyto water)^{-1}
% 
% % --------------------------------------------
% % LVH: (LVH hearts, Bache et al., Cardiovas Res, 1999)
% CT_Biopsy_LVH = 100.0; % umol (g dry wt)^{-1}
% ATPBL_Biopsy_LVH = 17.4; % umol (g dry wt)^{-1}
% CT_Biopsy_LVH_e = 3.2; % umol (g dry wt)^{-1}
% ATPBL_Biopsy_LVH_e = 0.5; % umol (g dry wt)^{-1}
% CrPATPBL_cellwater_LVH = 1.55; % Bache et al., Cardiovas Res, 1999
% CrPATPBL_cell_LVH = 1.55; % Bache et al., Cardiovas Res, 1999
% CrPATPBL_cell_LVH_e = 0.077; % Bache et al., Cardiovas Res, 1999
% x_AtC_LVH = 0.458e-3; % Corrected basal ATP consumption rate, mol s^{-1} (l cell)^{-1}
% 
% CT_cytowater_LVH = CT_Biopsy_LVH*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater_LVH = ATPBL_Biopsy_LVH*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% 
% CT_cellwater_LVH = CT_Biopsy_LVH*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater_LVH = ATPBL_Biopsy_LVH*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% 
% CT_cell_LVH = CT_Biopsy_LVH*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% ATPBL_cell_LVH = ATPBL_Biopsy_LVH*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% 
% CT_cytowater_LVH_e = CT_Biopsy_LVH_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater_LVH_e = ATPBL_Biopsy_LVH_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% CT_cellwater_LVH_e = CT_Biopsy_LVH_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater_LVH_e = ATPBL_Biopsy_LVH_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% CT_cell_LVH_e = CT_Biopsy_LVH_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% ATPBL_cell_LVH_e = ATPBL_Biopsy_LVH_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% 
% CrPBL_cell_LVH = CrPATPBL_cell_LVH*ATPBL_cell_LVH; % mol (l cell water)^{-1}
% CrPBL_cellwater_LVH = CrPATPBL_cellwater_LVH*ATPBL_cellwater_LVH; % mol (l cell water)^{-1}
% CrPBL_cytowater_LVH = CrPBL_cellwater_LVH*W_cell*rou_cell*Vcell_tissue/W_cyto/rou_cyto/Vcyto_tissue; % mol (l cyto water)^{-1}
% 
% % --------------------------------------------
% % Severe LVH: (LVH hearts, Zhang et al., JCI, 1993)
% CT_Biopsy_sLVH = 89.5; % umol (g dry wt)^{-1}
% ATPBL_Biopsy_sLVH = 14.7; % umol (g dry wt)^{-1}
% CT_Biopsy_sLVH_e = 5.0; % umol (g dry wt)^{-1}
% ATPBL_Biopsy_sLVH_e = 0.42; % umol (g dry wt)^{-1}
% CrPATPBL_cellwater_sLVH = 1.51; % Zhang et al., JCI, 1993
% CrPATPBL_cell_sLVH = 1.51; % Zhang et al., JCI, 1993
% CrPATPBL_cell_sLVH_e = 0.11; % Zhang et al., JCI, 1993
% x_AtC_sLVH = 0.458e-3; % Basal ATP consumption rate, mol s^{-1} (l cell)^{-1}, assumped to be same as Bache et al., Cardiovas Res, 1999
% 
% CT_cytowater_sLVH = CT_Biopsy_sLVH*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater_sLVH = ATPBL_Biopsy_sLVH*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% 
% CT_cellwater_sLVH = CT_Biopsy_sLVH*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater_sLVH = ATPBL_Biopsy_sLVH*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% 
% CT_cell_sLVH = CT_Biopsy_sLVH*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% ATPBL_cell_sLVH = ATPBL_Biopsy_sLVH*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% 
% CT_cytowater_sLVH_e = CT_Biopsy_sLVH_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% ATPBL_cytowater_sLVH_e = ATPBL_Biopsy_sLVH_e*(1-W_tissue)/Vcyto_tissue/rou_cyto/W_cyto*rou_water*1e-3; % mol (l cyto water)^{-1}
% CT_cellwater_sLVH_e = CT_Biopsy_sLVH_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% ATPBL_cellwater_sLVH_e = ATPBL_Biopsy_sLVH_e*(1-W_tissue)/Vcell_tissue/rou_cell/W_cell*rou_water*1e-3; % mol (l cell water)^{-1}
% CT_cell_sLVH_e = CT_Biopsy_sLVH_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% ATPBL_cell_sLVH_e = ATPBL_Biopsy_sLVH_e*(1-W_tissue)/Vcell_tissue*1e-3; % mol (l cell water)^{-1}
% 
% CrPBL_cell_sLVH = CrPATPBL_cell_sLVH*ATPBL_cell_sLVH; % mol (l cell water)^{-1}
% CrPBL_cellwater_sLVH = CrPATPBL_cellwater_sLVH*ATPBL_cellwater_sLVH; % mol (l cell water)^{-1}
% CrPBL_cytowater_sLVH = CrPBL_cellwater_sLVH*W_cell*rou_cell*Vcell_tissue/W_cyto/rou_cyto/Vcyto_tissue; % mol (l cyto water)^{-1}
% 
% %% For the single heart data points (Zhang et al., J Clin. Invest., 1993)
% % -------------------------------------------------------------
% % Set up CRtot and basal ATP for different experimental data points 
% % Assumptions of linear relationship between basal ATP and LVW/BW, total CrP and LVW/BW
% 
% % Loading the data
% LVWBW_con = 4.3; LVWBW_con_e = 0.14;
% LVWBW_mLVH = 5.7; LVWBW_mLVH_e = 0.19;
% LVWBW_LVH = 9.35; LVWBW_LVH_e = 0.97;
% LVWBW_sLVH = 8.9; LVWBW_sLVH_e = 0.67;
% 
% % Individual hearts: Taken from Figure 6, Zhang et al., J Clin. Invest., 1993
% CrPATPBL_exp_LVH_single = [1.99 1.93 1.8 1.54 1.44 1.22]; % single LVH heart
% LVWBW_exp_LVH_single =  [5.01 6.01 7.04 8.21 9.01 9.21]; % single LVH heart
% CrPATPBL_exp_LVHDD_single = [1.77 1.34 1.14 1.61 1.22 1.15]; % single LVH heart with diastolic dysfunction (filling pressure > 15 mmHg)
% LVWBW_exp_LVHDD_single =  [7.6 8.75 10.7 11.5 11.9 12.2]; % single LVH heart with diastolic dysfunction (filling pressure > 15 mmHg)
% % Combine the data from the single hearts
% CrPATPBL_exp_single = [CrPATPBL_exp_LVH_single CrPATPBL_exp_LVHDD_single];
% LVWBW_exp_single = [LVWBW_exp_LVH_single LVWBW_exp_LVHDD_single];
% 
% % Grouped hearts: Taken from Zhang et al.'s previous work
% CrPATPBL_exp_group = [CrPATPBL_cell CrPATPBL_cell_mLVH CrPATPBL_cell_LVH CrPATPBL_cell_sLVH]; % average values
% CrPATPBL_exp_group_e = [CrPATPBL_cell_e CrPATPBL_cell_mLVH_e CrPATPBL_cell_LVH_e CrPATPBL_cell_sLVH_e]; % standard deviations
% LVWBW_exp_group  =  [LVWBW_con LVWBW_mLVH LVWBW_LVH LVWBW_sLVH]; % average values
% LVWBW_exp_group_e =  [LVWBW_con_e  LVWBW_mLVH_e LVWBW_LVH_e LVWBW_sLVH_e]; % standard deviations
% CT_cytowater_group = [CT_cytowater CT_cytowater_mLVH CT_cytowater_LVH CT_cytowater_sLVH];
% CT_cellwater_group = [CT_cellwater CT_cellwater_mLVH CT_cellwater_LVH CT_cellwater_sLVH];
% CT_cellwater_group_e = [CT_cellwater_e CT_cellwater_mLVH_e CT_cellwater_LVH_e CT_cellwater_sLVH_e];
% CT_cell_group = CT_cellwater_group*rou_cell*W_cell/rou_water;
% CT_cell_group_e = CT_cellwater_group_e*rou_cell*W_cell/rou_water;
% ATPBL_cell_group = [ATPBL_cell ATPBL_cell_mLVH ATPBL_cell_LVH ATPBL_cell_sLVH];
% ATPBL_cell_group_e = [ATPBL_cell_e ATPBL_cell_mLVH_e ATPBL_cell_LVH_e ATPBL_cell_sLVH_e];
% CrPBL_cell_group = [CrPBL_cell CrPBL_cell_mLVH CrPBL_cell_LVH CrPBL_cell_sLVH];
% 
% % Compute CRtot, basal CrP and basal ATP (per cell) for the single hearts
% % Creatine pool per cytowater: assuming linear changes of CRtot versus LVW/BW
% % Find the linear relationship between CRtot and LVWBW based on data of
% % grouped heart
% pf_x1 = polyfit(LVWBW_exp_group([1 4]),CT_cytowater_group([1 4]),1);
% CT_cytowater_exp_single = pf_x1(1)*LVWBW_exp_single + pf_x1(2);
% CT_cellwater_exp_single = CT_cytowater_exp_single*Vcyto_tissue*rou_cyto*W_cyto/Vcell_tissue/rou_cell/W_cell;
% CT_cell_exp_single = CT_cellwater_exp_single*rou_cell*W_cell/rou_water;
% % CT_cytowater_exp_single1 = CT_cytowater_sLVH + (CT_cytowater-CT_cytowater_sLVH)*(LVWBW_exp_single-LVWBW_sLVH)/(LVWBW_con-LVWBW_sLVH);
% % CT_cellwater_exp_single1 = CT_cellwater_sLVH +
% % (CT_cellwater-CT_cellwater_sLVH)*(LVWBW_exp_single-LVWBW_sLVH)/(LVWBW_con-LVWBW_sLVH);
% % Basal ATP per cell: assuming linear changes of basal ATP versus LVW/BW
% % Find the linear relationship between ATPBL and LVWBW based on data of
% % grouped heart
% pf_x2 = polyfit(LVWBW_exp_group([1 4]),ATPBL_cell_group([1 4]),1);
% ATPBL_cell_exp_single = pf_x2(1)*LVWBW_exp_single + pf_x2(2);
% % ATPBL_cell_exp_single = ATPBL_cell_sLVH + (ATPBL_cell-ATPBL_cell_sLVH)*(LVWBW_exp_single-LVWBW_sLVH)/(LVWBW_con-LVWBW_sLVH); 
% % Basal CrP per cell: computed from basal CrP/ATP and basal ATP
% CrPBL_cell_exp_single = CrPATPBL_exp_single.*ATPBL_cell_exp_single;
% 
% % Combine the single and group hearts
% LVWBW_exp = [LVWBW_exp_single LVWBW_exp_group];
% CrPATPBL_exp = [CrPATPBL_exp_single CrPATPBL_exp_group];
% CT_cytowater_exp = [CT_cytowater_exp_single CT_cytowater_group];
% CT_cellwater_exp = [CT_cellwater_exp_single CT_cellwater_group];
% CT_cell_exp = CT_cellwater_exp*rou_cell*W_cell/rou_water;
% ATPBL_cell_exp = [ATPBL_cell_exp_single ATPBL_cell_group];
% CrPBL_cell_exp = [CrPBL_cell_exp_single CrPBL_cell_group];
% i_con = 13:16;
% i_LVH = 1:6;
% i_LVH_DD = 7:12;


%% Model Parameters
load param14.mat;     % Load model parameters of the previously-developed mitochondrial model (Wu et al., J Biol Chem 282 (34), 24525 (2007)
param(39) = 50;       % Decrease leaking activity of proton for in vivo mito
param(35) = 6.762e-3; % Reset ANT activity for an updated ANT kinetic module 
parma(42) = 1e-3;     % Ref_PYR
param(43) = 4e-4;     % Carnitine in the Cytosol (400 uM)
