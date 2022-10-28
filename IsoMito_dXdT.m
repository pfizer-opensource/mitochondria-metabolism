function [f] = IsoMito_dXdT(t,x,param,ExpType,StateType,J_AtC)

% Computing the time derivatives of state variables, according to types of
% experiments or energy states

% Implementing computer model presented in Wu et al.'s "Computer modeling of
% mitochondrial TCA cycle, oxidative phosphorylation, metabolite transport, and
% electrophysiology" and its Appendices

% Update from Version: 1.4 (previous update on 10/27/08)

%% Defining global variables
global iH_x idPsi iATP_x iADP_x iAMP_x iGTP_x iGDP_x iPi_x iNADH_x ...
       iQH2_x iOAA_x iACCOA_x iCIT_x iICIT_x iAKG_x iSCOA_x iCOASH_x ...
       iSUC_x iFUM_x iMAL_x iGLU_x iASP_x iK_x iMg_x iO2_x iCO2tot_x ...
       iCred_i iATP_i iADP_i iAMP_i iPi_i iH_i iMg_i iK_i ...
       iATP_c iADP_c iPi_c iH_c iMg_c iK_c ...
       iPYR_x iPYR_i iPYR_c iCIT_i iCIT_c iAKG_i iAKG_c iSUC_i iSUC_c ...
       iMAL_i iMAL_c iASP_i iASP_c iGLU_i iGLU_c iFUM_i iFUM_c ...
       iICIT_i iICIT_c iGLC_c iG6P_c iPCr_c iAMP_c;
     
%% Defining indices for all reactions in the TCA cycle
ipdh    = 1;
icits   = 2;
iacon   = 3;
iisod   = 4;
iakgd   = 5;
iscoas  = 6;
isdh    = 7;
ifum    = 8;
imdh    = 9;
indk    = 10;
igot    = 11;
%% Setting adjustable parameter values
% Vmax values of TCA cycle fluxes
Vmax(1) = param(1);
Vmax(2) = param(2);
Vmax(3) = param(3);
Vmax(4) = param(4);
Vmax(5) = param(5);
Vmax(6) = param(6);
Vmax(7) = param(7);
Vmax(8) = param(8);
Vmax(9) = param(9);
Vmax(10)= param(10);
Vmax(11)= param(11);

% Activities of transporters
x_PYR_H   = param(12);
x_GLU_H   = param(13);
x_CIT_MAL = param(14); 
x_AKG_MAL = param(15); 
% param(16) not used
x_MAL_PI  = param(17); 
x_ASP_GLU = param(18);
% param(19) not used
% param(20) not used
x_SUC_MAL = param(21); 
% param(22) not used
% param(23) was used by Kir1 in J_akgd
% param(24) was not used 

% Hexokinase activity (for LaNoue state-3 experiments)
if (ExpType == 1) && (StateType == 3)
    x_HK = param(25);        
else 
    x_HK = 0;
end


% Parameters for oxidative phosphorylation
x_C1 = param(31);
x_C3 = param(32);
x_C4 = param(33);
x_F1 = param(34);
x_ANT = param(35);
x_Pi1 = param(36);
k_PiH = param(37);
x_KH = param(38);
x_Hle = param(39);
k_Pi1 = param(40);
k_Pi2 = param(41);

% Permeability coefficients for TCA intermediates
p_TI     = 85; % assumed to be equal to x_A for nucleotides (micro sec^-1)
x_PYRt   = p_TI;
x_GLUt   = p_TI;
x_ASPt   = p_TI;
x_CITt   = p_TI;
x_ICITt  = p_TI;
x_AKGt   = p_TI;
x_FUMt   = 0*p_TI;
x_SUCt   = p_TI;
x_MALt   = p_TI;

   
%% Defining indices of reactants
iH      =   1;
iATP    =   2;
iADP    =   3;
iAMP    =   4;
iGTP    =   5;
iGDP    =   6;
iPi     =   7;
iNADH   =   8;
iNAD    =   9;
iQH2    =   10;
iCOQ    =   11;
iOAA    =   12;
iACCOA  =   13;
iCIT    =   14;
iICIT   =   15;
iAKG    =   16;
iSCOA   =   17;
iCOASH  =   18;
iSUC    =   19;
iFUM    =   20;
iMAL    =   21;
iGLU    =   22;
iASP    =   23;
iK      =   24;
iMg     =   25;
iCox    =   26;
iCred   =   27;
iO2     =   28;
iH2O    =   29;
iFADH2  =   30;
iFAD    =   31;
iCO2tot =   32;
iPCr    =   33;
iCr     =   34;
iPYR    =   35;
iGLC    =   36;
iG6P    =   37;

N_reactant = 37;

%%  Listing fixed model parameters
%  (i) Thermochemical constants
if ExpType == 1
    RT = 8.314*(28+273.15)/1e3;          % kJ  mol^{-1}
else
    RT = 8.314*(37+273.15)/1e3;          % kJ  mol^{-1}
end;
F = 0.096484;                   		 % kJ mol^{-1} mV^{-1}

%  (ii) Subcelular volumes and water spaces
if (ExpType == 1) % LaNoue et al.'s experiments (JBC, 1970)
    W_c = 80; % buffer water space for LaNoue experiments
elseif (ExpType == 2)  % Bose et al.'s experiments (JBC, 2003)
    W_c = 80*10; % buffer water space for Bose experiments
elseif (ExpType == 3)|| (ExpType == 8)   % Caridac Tissue
    Vmito = 0.2882;                 % (ml mito / ml cell)
    Vcyto = 0.6801;                 % (ml cyto / ml cell)
    Rm_cyto = Vmito / Vcyto;        % Volume ratio mito volume / cytoplasm volume
    Rm_cell = Vmito;                % Volume ratio mito volume / cell volume
    Rc_cell = Vcyto;                % Volume ratio of cytosol / cell volume
    W_c = 0.807*1.044;              % cytosol water space (ml water per ml cytosol) [VB2002]
elseif (ExpType == 4) % Skeletal Muscle
    Vmito = 0.056;                  % (ml mito / ml cell)
    Vcyto = 0.894;                  % (ml cyto / ml cell)
    Rm_cyto = Vmito / Vcyto;        % Volume ratio mito volume / cytoplasm volume
    Rm_cell = Vmito;                % Volume ratio mito volume / cell volume
    Rc_cell = Vcyto;                % Volume ratio of cytosol / cell volume
    W_c = 0.807*1.044;              % cytosol water space (ml water per ml cytosol) [VB2002]
elseif (ExpType == 7)
   rho_m = 3.6697e-6;  % (l mito) (mg protein)^{-1} %%assume 0.025 mg/mL * 2mL rxn volume
    W_c = 1090; % 2 mL *0.025 mg protein ref.supp. Muoio 
else
    W_c = 80;
end
W_m = 0.664*1.09;               % mitochondrial water space (ml water per ml mito) [VB2002]
W_x = 0.9*W_m;                  % Matrix water space  (90% is matrix)
W_i = 0.1*W_m;                  % IM water space      (10% is intermembranee-space)

%  (iii) Pooled concentrations
Ctot   = 2.70e-3;               % M; total cytoC
Qtot   = 1.35e-3;               % M; total Q + QH2
NADtot = 2.97e-3;               % M; total NAD+NADH
FADtot = 0.1e-3;                % M; total FADred+FADox (from Feng Yang)

%  (iv) Ox Phos Model Parameters
n_A     = 3.0;   % unitless  3 protons
k_O2   = 1.2e-4; % kinetic constant for complex 4
k_mADP = 3.5e-6; % ANT MM constant

%  (v) Outer membrane transport parameters
x_A    = 85;        % micron sec^{-1}
x_Pi2  = 327;       % micron sec^{-1}
gamma  = 5.99;      % mito membrane area per cell volume micron^{-1}



%% Loading values of state variables
MinCon = 1e-32;

% (i) Matrix species and dPsi
dPsi        = x(idPsi);
ATP_x       = x(iATP_x);
ADP_x       = x(iADP_x);
AMP_x       = x(iAMP_x);
GTP_x       = x(iGTP_x);
GDP_x       = x(iGDP_x);
Pi_x        = x(iPi_x);
NADH_x      = x(iNADH_x);
QH2_x       = x(iQH2_x);
PYR_x       = x(iPYR_x);
OAA_x       = x(iOAA_x);
ACCOA_x     = x(iACCOA_x);
CIT_x       = x(iCIT_x);
ICIT_x      = x(iICIT_x);
AKG_x       = x(iAKG_x);
SCOA_x      = x(iSCOA_x);
COASH_x     = x(iCOASH_x);
SUC_x       = x(iSUC_x);
FUM_x       = x(iFUM_x);
MAL_x       = x(iMAL_x);
GLU_x       = x(iGLU_x);
ASP_x       = x(iASP_x);
H_x         = x(iH_x); 
K_x         = x(iK_x);
Mg_x        = x(iMg_x);
O2          = x(iO2_x);
CO2tot      = x(iCO2tot_x);

% (ii) IM space species
Cred_i      = max(0,x(iCred_i));
ATP_i       = x(iATP_i);
ADP_i       = x(iADP_i);
AMP_i       = x(iAMP_i);
Pi_i        = x(iPi_i);
PYR_i       = x(iPYR_i);
CIT_i       = x(iCIT_i);
AKG_i       = x(iAKG_i);
SUC_i       = x(iSUC_i);
MAL_i       = x(iMAL_i);
GLU_i       = x(iGLU_i);
ASP_i       = x(iASP_i);
H_i         = x(iH_i); 
Mg_i        = x(iMg_i);
K_i         = x(iK_i);
FUM_i       = x(iFUM_i);
ICIT_i      = x(iICIT_i);

% (iii) Cytoplasmic species
ATP_c       = x(iATP_c);
ADP_c       = x(iADP_c);
Pi_c        = x(iPi_c);
PYR_c       = x(iPYR_c);
CIT_c       = x(iCIT_c);
AKG_c       = x(iAKG_c);
SUC_c       = x(iSUC_c);
MAL_c       = x(iMAL_c);
GLU_c       = x(iGLU_c);
ASP_c       = x(iASP_c);
H_c         = x(iH_c);
Mg_c        = x(iMg_c);
K_c         = x(iK_c);
FUM_c       = x(iFUM_c);
ICIT_c      = x(iICIT_c);
GLC_c       = x(iGLC_c);
G6P_c       = x(iG6P_c);
AMP_c       = x(iAMP_c);

% (iv) Other concentrations computed from the state variables:
NAD_x  = NADtot - NADH_x;
COQ_x  = Qtot - QH2_x;
Cox_i  = Ctot - Cred_i;

% (v) set the H+, Mg2+, and K+ are permeable for outer mito membrane
H_i = H_c;
Mg_i = Mg_c;
K_i = K_c;


%% Loading thermodynamic data (deltG, pK, etc.)
% T = 298.15K (25 C)    I = 0.17 M
% standard Gibbs free energy of formation of reference species, (kJ/mol)

% without temperature correction on dGf
dGf1(1:N_reactant) = 0;
dGf1(iH2O) = -235.74;        % H2O
dGf1(iO2) = 16.40;           % O2(aq)
dGf1(iNADH) = 39.31;         % NADH
dGf1(iNAD) = 18.10;          % NAD+
dGf1(iQH2) = -23.30;         % QH2
dGf1(iCOQ) = 65.17;          % Q
dGf1(iATP) = -2771.00;       % ATP4-
dGf1(iADP) = -1903.96;       % ADP3-
dGf1(iAMP) = -1034.66;       % AMP2-
dGf1(iGTP) = dGf1(iATP);
dGf1(iGDP) = dGf1(iADP);
dGf1(iCred) = -27.41;        % CytoC(red)2+
dGf1(iCox) = -6.52;          % CytoC(ox)3+
dGf1(iPi) = -1098.27;        % HPO42-
dGf1(iPCr) = 0;              % PCr2-
dGf1(iCr) = -252.68;         % HCr
dGf1(iFADH2) = -67.60;       % FADH2-enz
dGf1(iFAD) = 19.55;          % FAD-enz
dGf1(iCOASH) = -0.72;        % CoAS-
dGf1(iACCOA) = -178.19;      % AcCoA
dGf1(iOAA) = -794.74;        % OAA2-
dGf1(iCIT) = -1165.59;       % CIT3-
dGf1(iICIT) = -1158.94;      % ICIT3-
dGf1(iAKG) = -793.41;        % AKG2-
dGf1(iSCOA) = -507.55;       % SCoA-
dGf1(iSUC) = -690.44;        % SUC2-
dGf1(iFUM) = -603.32;        % FUM2-
dGf1(iMAL) = -842.66;        % MAL2-
dGf1(iASP) = -692.26;        % ASP-
dGf1(iGLU) = -692.40;        % GLU- (L-glutamate)
dGf1(iCO2tot) = -530.71;     % CO2tot
dGf1(iPYR) = -470.82;        % PYR2-
dGf1(iGLC) = -907.21;        % Glucose
dGf1(iG6P) = -1758.87;       % Glucose-6-phosphate

% K values for reference species 
% pK_KATP is corrected to be 1.013, 08/26/08
% pK_KADP is corrected to be 0.882, 08/26/08
% pK_KAMP is corrected to be 0.6215, 08/26/08
% pK_MgOAA is corrected to be 0.8629, 08/26/08
% pK_KSUC is corrected to be 0.3525, 08/26/08
Kh(1:N_reactant) = inf; Km(1:N_reactant) = inf; Kk(1:N_reactant) = inf;
Kh(iATP) = 10^(-6.59); Km(iATP) = 10^(-3.82); Kk(iATP) = 10^(-1.013); 
Kh(iADP) = 10^(-6.42); Km(iADP) = 10^(-2.79); Kk(iADP) = 10^(-0.882); 
Kh(iAMP) = 10^(-6.22); Km(iAMP) = 10^(-1.86); Kk(iAMP) = 10^(-0.6215); 
Kh(iGTP) = Kh(iATP); Km(iGTP) = Km(iATP); Kk(iGTP) = Kk(iATP);
Kh(iGDP) = Kh(iADP); Km(iGDP) = Km(iADP); Kk(iGDP) = Kk(iADP);
Kh(iPi) = 10^(-6.71); Km(iPi) = 10^(-1.69); Kk(iPi) = 10^(+0.0074);
Kh(iCOASH) = 10^(-8.13);
Km(iOAA) = 10^(-0.8629);                                                
Kh(iCIT) = 10^(-5.63); Km(iCIT) = 10^(-3.37); Kk(iCIT) = 10^(-0.339);
Kh(iICIT) = 10^(-5.64); Km(iICIT) = 10^(-2.46); 
Kh(iSCOA) = 10^(-3.96);
Kh(iSUC) = 10^(-5.13); Km(iSUC) = 10^(-1.17); Kk(iSUC) = 10^(-0.3525);   
Kh(iFUM) = 10^(-4.10);
Kh(iMAL) = 10^(-4.75); Km(iMAL) = 10^(-1.55); Kk(iMAL) = 10^(+0.107);
Kh(iCO2tot) = 10^(-9.82);
Km(iPYR) = 10^(-1.02); 
Kh(iG6P) = 10^(-5.91);
% Kh(iGLU) = 10^(-4.25);   % from Nelson & Cox, "Lehninger's Princinples of Biochemistry", p78      
Kh(iGLU) = 10^(-4.06);  % 37 C, I = 0.15
Km(iGLU) = 10^(-1.82);
% Kh(iASP) = 10^(-3.65);   % from Nelson & Cox, "Lehninger's Princinples of Biochemistry", p78 
Kh(iASP) = 10^(-3.65); % 37 C, I = 0.15
Km(iASP) = 10^(-2.32);

% compute binding polynomials for reactants
P_x(1:N_reactant) = 1; P_c(1:N_reactant) = 1; P_i(1:N_reactant) =1;
P_x(iATP) = 1 + H_x/Kh(iATP) + Mg_x/Km(iATP) + K_x/Kk(iATP);
P_c(iATP) = 1 + H_c/Kh(iATP) + Mg_c/Km(iATP) + K_c/Kk(iATP);
P_i(iATP) = 1 + H_i/Kh(iATP) + Mg_i/Km(iATP) + K_i/Kk(iATP);
P_x(iADP) = 1 + H_x/Kh(iADP) + Mg_x/Km(iADP) + K_x/Kk(iADP);
P_c(iADP) = 1 + H_c/Kh(iADP) + Mg_c/Km(iADP) + K_c/Kk(iADP);
P_i(iADP) = 1 + H_i/Kh(iADP) + Mg_i/Km(iADP) + K_i/Kk(iADP);
P_x(iAMP) = 1 + H_x/Kh(iAMP) + Mg_x/Km(iAMP) + K_x/Kk(iAMP);
P_c(iAMP) = 1 + H_c/Kh(iAMP) + Mg_c/Km(iAMP) + K_c/Kk(iAMP);
P_i(iAMP) = 1 + H_i/Kh(iAMP) + Mg_i/Km(iAMP) + K_i/Kk(iAMP);
P_x(iGTP) = P_x(iATP); 
P_c(iGTP) = P_c(iATP);
P_i(iGTP) = P_i(iATP);
P_x(iGDP) = P_x(iADP); 
P_c(iGDP) = P_c(iADP);
P_i(iGDP) = P_i(iADP);
P_x(iPi) = 1 + H_x/Kh(iPi) + Mg_x/Km(iPi) + K_x/Kk(iPi); % add K-bound item, 06/10/08
P_c(iPi) = 1 + H_c/Kh(iPi) + Mg_c/Km(iPi) + K_c/Kk(iPi); % add K-bound item, 06/10/08
P_i(iPi) = 1 + H_i/Kh(iPi) + Mg_i/Km(iPi) + K_i/Kk(iPi); % add K-bound item, 06/10/08
P_x(iCOASH) = 1 + H_x/Kh(iCOASH);
P_i(iCOASH) = 1 + H_i/Kh(iCOASH);
P_c(iCOASH) = 1 + H_c/Kh(iCOASH);
P_x(iOAA) = 1 + Mg_x/Km(iOAA);
P_i(iOAA) = 1 + Mg_i/Km(iOAA);
P_c(iOAA) = 1 + Mg_c/Km(iOAA);
P_x(iCIT) = 1 + H_x/Kh(iCIT) + Mg_x/Km(iCIT) + K_x/Kk(iCIT);
P_i(iCIT) = 1 + H_i/Kh(iCIT) + Mg_i/Km(iCIT) + K_i/Kk(iCIT);
P_c(iCIT) = 1 + H_c/Kh(iCIT) + Mg_c/Km(iCIT) + K_c/Kk(iCIT);
P_x(iICIT) = 1 + H_x/Kh(iICIT) + Mg_x/Km(iICIT);
P_i(iICIT) = 1 + H_i/Kh(iICIT) + Mg_i/Km(iICIT);
P_c(iICIT) = 1 + H_c/Kh(iICIT) + Mg_c/Km(iICIT);
P_x(iSCOA) = 1 + H_x/Kh(iSCOA);
P_i(iSCOA) = 1 + H_i/Kh(iSCOA);
P_c(iSCOA) = 1 + H_c/Kh(iSCOA);
P_x(iSUC) = 1 + H_x/Kh(iSUC) + Mg_x/Km(iSUC) + K_x/Kk(iSUC);
P_i(iSUC) = 1 + H_i/Kh(iSUC) + Mg_i/Km(iSUC) + K_i/Kk(iSUC);
P_c(iSUC) = 1 + H_c/Kh(iSUC) + Mg_c/Km(iSUC) + K_c/Kk(iSUC);
P_x(iFUM) = 1 + H_x/Kh(iFUM);
P_i(iFUM) = 1 + H_i/Kh(iFUM);
P_c(iFUM) = 1 + H_c/Kh(iFUM);
P_x(iMAL) = 1 + H_x/Kh(iMAL) + Mg_x/Km(iMAL) + K_x/Kk(iMAL);
P_i(iMAL) = 1 + H_i/Kh(iMAL) + Mg_i/Km(iMAL) + K_i/Kk(iMAL);
P_c(iMAL) = 1 + H_c/Kh(iMAL) + Mg_c/Km(iMAL) + K_c/Kk(iMAL);
P_x(iCO2tot) = 1 + H_x/Kh(iCO2tot);
P_i(iCO2tot) = 1 + H_i/Kh(iCO2tot);
P_c(iCO2tot) = 1 + H_c/Kh(iCO2tot);
P_x(iPYR) = 1 + Mg_x/Km(iPYR);
P_i(iPYR) = 1 + Mg_i/Km(iPYR);
P_c(iPYR) = 1 + Mg_c/Km(iPYR);
P_x(iG6P) = 1 + H_x/Kh(iG6P);
P_i(iG6P) = 1 + H_i/Kh(iG6P);
P_c(iG6P) = 1 + H_c/Kh(iG6P);
P_x(iGLU) = 1 + H_x/Kh(iGLU) + Mg_x/Km(iGLU); % correct Mg-bound item, 06/10/08
P_i(iGLU) = 1 + H_i/Kh(iGLU) + Mg_i/Km(iGLU); % correct Mg-bound item, 06/10/08
P_c(iGLU) = 1 + H_c/Kh(iGLU) + Mg_c/Km(iGLU); % correct Mg-bound item, 06/10/08 
P_x(iASP) = 1 + H_x/Kh(iASP) + Mg_x/Km(iASP); % correct Mg-bound item, 06/10/08
P_i(iASP) = 1 + H_i/Kh(iASP) + Mg_i/Km(iASP); % correct Mg-bound item, 06/10/08
P_c(iASP) = 1 + H_c/Kh(iASP) + Mg_c/Km(iASP); % correct Mg-bound item, 06/10/08


%% I. Flux expresssions in the TCA cycle 

%% -------------------------------
% 1. Pyruvate dehydrogenase
% PYR + COASH + NAD (+H2O) = CO2tot + SCOA + NADH
% A - PYR; B - COASH; C - NAD; P - CO2tot; Q - SCOA; R - NADH;

% load concentrations of reactants and products
A = PYR_x;
B = COASH_x;
C = NAD_x;
P = CO2tot;
Q = ACCOA_x;
R = NADH_x;

% dG and Keq vlaues
dGr_pdho = dGf1(iCO2tot) + dGf1(iACCOA) + dGf1(iNADH) ...
          - dGf1(iPYR) - dGf1(iCOASH) - dGf1(iNAD) - dGf1(iH2O);
Keq_pdho = exp(-dGr_pdho/RT);
Keq_pdh = Keq_pdho*(1/H_x)*(P_x(iCO2tot)*P_x(iACCOA)*P_x(iNADH)) ...
                           /(P_x(iPYR)*P_x(iCOASH)*P_x(iNAD));

% Km and Ki values (Molar)
KmA = 38.3e-6; 
KmB = 9.9e-6;
KmC = 60.7e-6;
KiACCOA = 40.2e-6;
KiNADH = 40.0e-6;

% Inhibition Constants
ai1 = 1 + ACCOA_x/KiACCOA;
ai2 = 1 + NADH_x/KiNADH;

% Vm values
Vmf = Vmax(ipdh);                  

% total reaction flux
if (A > MinCon) && (B > MinCon) && (C > MinCon) 
    J_pdh = Vmf*(A*B*C-P*Q*R/Keq_pdh)/ (KmC*ai2*A*B + KmB*ai1*A*C + KmA*B*C + A*B*C);
else 
    J_pdh = 0;
end

%% -------------------------------
% 2. Citrate synthetase
% OAA + ACCOA (+H2O) = COASH + CIT
% A - OAA; B - ACCOA; P - COASH; Q - CIT;

% load concentrations of reactants and products
A = OAA_x;
B = ACCOA_x;
P = COASH_x;
Q = CIT_x;
    
% dG and Keq values
dGr_citso = dGf1(iCOASH) + dGf1(iCIT) - dGf1(iACCOA) - dGf1(iOAA) - dGf1(iH2O);
Keq_citso = exp(-dGr_citso/RT);
Keq_cits = Keq_citso*(1/H_x^2)*(P_x(iCOASH)*P_x(iCIT)) ...
                           /(P_x(iACCOA)*P_x(iOAA));
 
% Km and Ki values (Molar)
KmA = 4e-6; 
KmB = 14e-6;
Kia = 3.33e-6;
KiCIT = 1600e-6;
KiATP = 900e-6;
KiADP = 1800e-6;
KiAMP = 6000e-6;
KiCOASH = 67e-6;
KiSCOA = 140e-6;

% inhibition coefficients
uCIT_x = CIT_x * (1+H_x/Kh(iCIT))/P_x(iCIT); % unchelated
uATP_x = ATP_x * (1+H_x/Kh(iATP))/P_x(iATP); % unchelated
uADP_x = ADP_x * (1+H_x/Kh(iADP))/P_x(iADP); % unchelated
uAMP_x = AMP_x * (1+H_x/Kh(iAMP))/P_x(iAMP); % unchelated
ai1 = 1 + uCIT_x/KiCIT;
ai2 = 1 + uATP_x/KiATP + uADP_x/KiADP + uAMP_x/KiAMP ...
        + COASH_x/KiCOASH + SCOA_x/KiSCOA;
                      
% Vm values
Vmf = Vmax(icits);

% forward reaction flux
J_cits_f = Vmf*A*B / (Kia*KmB*ai1 + KmA*ai1*B + KmB*ai2*A + A*B);

% overall reaction flux
J_cits = J_cits_f - Vmf*(P*Q/Keq_cits) / (Kia*KmB*ai1 + KmA*ai1*B + KmB*ai2*A + A*B);


%% -------------------------------
% 3. Aconitase
% CIT = ICIT
% A - CIT; P - ICIT;

% load concentrations of reactants and products
A = CIT_x;
P = ICIT_x;
    
% dG and Keq values
dGr_acono = dGf1(iICIT) - dGf1(iCIT);
Keq_acono = exp(-dGr_acono/RT);
% Keq_acon = Keq_acono;
Keq_acon = Keq_acono*P_x(iICIT)/P_x(iCIT);
           
% Km and Ki values (Molar)
KmA = 1161e-6; 
KmP = 434e-6;
            
% Vm values
Vmf = Vmax(iacon);
Vmr = Vmf*(KmP/KmA/Keq_acon);

% forward reaction flux
J_acon_f = Vmf*Vmr*A /(KmA*Vmr+Vmr*A+Vmf/Keq_acon*P);

% total reaction flux
J_acon = J_acon_f - Vmf*Vmr*(P/Keq_acon)/(KmA*Vmr+Vmr*A+Vmf/Keq_acon*P);


%% -------------------------------
% 4. Isocitrate dehydrogenase
% NAD + ICIT (+ H2O) =  AKG + NADH + CO2tot
% A - NAD; B - ICIT; P - AKG; Q - NADH; R - CO2tot;

% load concentrations of reactants and products
A = NAD_x;
B = ICIT_x;
P = AKG_x;
Q = NADH_x;
R = CO2tot;
    
% dG and Keq values
dGr_isodo = dGf1(iAKG) + dGf1(iNADH) + dGf1(iCO2tot) ...
            - dGf1(iICIT) - dGf1(iNAD) - dGf1(iH2O);
Keq_isodo = exp(-dGr_isodo/RT);
% Keq_isod = Keq_isodo*(1/H_x^2);
Keq_isod = Keq_isodo*(1/H_x^2)*(P_x(iAKG)*P_x(iNADH)*P_x(iCO2tot)) ...
                           /(P_x(iICIT)*P_x(iNAD));
                       
% Km and Ki values (Molar)
KmA = 74e-6;
KmB = 183e-6;
nH = 3.0;
Kib = 23.8e-6;
Kiq = 29e-6;
KiATP = 91e-6;
KaADP = 50e-6;

% inhibition coefficients
fATP_x = ATP_x * (1+H_x/Kh(iATP))/P_x(iATP);
fADP_x = ADP_x * (1+H_x/Kh(iADP))/P_x(iADP);
ai = 1 + KaADP/fADP_x*(1+fATP_x/KiATP);

% Vm values
Vmf = Vmax(iisod);

% total reaction flux
if (A > MinCon) && (B > MinCon) 
    J_isod = Vmf/(1+(KmB/B)^nH*ai+KmA/A*(1+(Kib/B)^nH*ai+Q*ai/Kiq))*(1-1/Keq_isod*P*Q*R/A/B);
else
    J_isod = 0;
end


%% -------------------------------
% 5. alpha-Ketoglutarate dehydrogenase
% AKG + COASH + NAD (+ H2O) = CO2tot + SCOA + NADH
% A - AKG; B - COASH; C - NAD; P - CO2tot; Q - SCOA; R - NADH;

% load concentrations of reactants and products
A = AKG_x;
B = COASH_x;
C = NAD_x;
P = CO2tot;
Q = SCOA_x;
R = NADH_x;
    
% dG and Keq values
dGr_akgdo = dGf1(iCO2tot) + dGf1(iSCOA) + dGf1(iNADH) ...
            - dGf1(iAKG) - dGf1(iCOASH) - dGf1(iNAD) - dGf1(iH2O);
Keq_akgdo = exp(-dGr_akgdo/RT);
% Keq_akgd = Keq_akgdo*(1/H_x);
Keq_akgd = Keq_akgdo*(1/H_x)*(P_x(iCO2tot)*P_x(iSCOA)*P_x(iNADH)) ...
                           /(P_x(iAKG)*P_x(iCOASH)*P_x(iNAD));
                       
% Km and Ki values (Molar)
KmA = 80e-6;
KmB = 55e-6;
KmC = 21e-6;
Kiq = 6.9e-6;
Kir1 = 4.5e-6;
Kir2 = 12.7e-6; 
KiATP = 50e-6;
KaADP = 100e-6;

% inhibition coefficients
fATP_x = ATP_x * (1+H_x/Kh(iATP))/P_x(iATP);
fADP_x = ADP_x * (1+H_x/Kh(iADP))/P_x(iADP);
ai = 1 + KaADP/fADP_x*(1+fATP_x/KiATP);
                       
% Vm values
Vmf = Vmax(iakgd);

% % for testing
Kir1 = param(23);
Kir2 = 1e3;

% total reaction flux
if (A > MinCon) && (B > MinCon) && (C > MinCon) 
    J_akgd = Vmf/(1+KmA/A*ai+KmB/B*(1+Q/Kiq)+KmC/C*(1+R/Kir1))/(1+R/Kir2)*(1-1/Keq_akgd*P*Q*R/A/B/C);
else
    J_akgd = 0;
end

 
%% -------------------------------
% 6. Succinyl-CoA synthetase
% GDP + SCOA + PI = COASH + SUC + GTP

% load concentrations of reactants and products
A = GDP_x;
B = SCOA_x;
C = Pi_x;
P = COASH_x;
Q = SUC_x;
R = GTP_x;

    
% dG and Keq values
dGr_scoaso = dGf1(iCOASH) + dGf1(iSUC) + dGf1(iGTP) ...
            - dGf1(iGDP) - dGf1(iSCOA) - dGf1(iPi);
Keq_scoaso = exp(-dGr_scoaso/RT);
% Keq_scoas = Keq_scoaso*(1/H_x);
Keq_scoas = Keq_scoaso*(1/H_x)*(P_x(iCOASH)*P_x(iSUC)*P_x(iGTP)) ...
                           /(P_x(iSCOA)*P_x(iPi)*P_x(iGDP));

% Km and Ki values (Molar)
KmA = 16e-6;
KmB = 55e-6;
KmC = 660e-6;
KmP = 20e-6;
KmQ = 880e-6;
KmR = 11.1e-6;
Kia = 5.5e-6;
Kib = 100e-6;
Kic = 2000e-6;
Kip = 20e-6;
Kiq = 3000e-6;
Kir = 11.1e-6;
                       
% Vm values
Vmf = Vmax(iscoas);
Vmr = Vmf/Keq_scoas*KmP*Kiq*Kir/(Kia*Kib*KmC);

% total reaction flux  
% Correct typos in J_scoas expression (10/27/08): 
% Vmf*KmQ*Kir*A*B*Q ->Vmf*KmQ*Kir*A*B*P; 
% Vmf*KmA*B*C*Q*R -> Vmr*KmA*B*C*Q*R;
% Vmf*KmA*B*C*P*Q*R -> Vmr*KmA*B*C*P*Q*R
J_scoas = (Vmf*Vmr*A*B*C - Vmf*Vmr*(P*Q*R/Keq_scoas)) ...
  /(Vmr*Kia*Kib*KmC+Vmr*Kib*KmC*A+Vmr*Kia*KmB*C ...
    + Vmr*KmC*A*B+Vmr*KmB*A*C+Vmr*KmA*B*C+Vmr*A*B*C ...
    + Vmf*Kir*KmQ*P/Keq_scoas+Vmf*Kiq*KmP*R/Keq_scoas+Vmf*KmR*P*Q/Keq_scoas+Vmf*KmQ*P*R/Keq_scoas...
    + Vmf*KmP*Q*R/Keq_scoas+Vmf*P*Q*R/Keq_scoas+Vmf*KmQ*Kir*A*P/Kia/Keq_scoas+Vmr*Kia*KmB*C*R/Kir...
    + Vmf*KmQ*Kir*A*B*P/Kia/Kib/Keq_scoas+Vmr*KmA*B*C*R/Kir+Vmf*KmR*A*P*Q/Kia/Keq_scoas...
    + Vmr*Kia*KmB*C*Q*R/Kiq/Kir+Vmf*Kir*KmQ*A*B*C*P/Kia/Kib/Kic/Keq_scoas+Vmf*Kip*KmR*A*B*C*Q/Kia/Kib/Kic/Keq_scoas...
    + Vmf*KmR*A*B*P*Q/Kia/Kib/Keq_scoas+Vmr*KmA*B*C*Q*R/Kiq/Kir+Vmr*KmA*Kic*B*P*Q*R/Kip/Kiq/Kir...
    + Vmr*Kia*KmB*C*P*Q*R/Kip/Kiq/Kir+Vmf*KmR*A*B*C*P*Q/Kia/Kib/Kic/Keq_scoas+Vmr*KmA*B*C*P*Q*R/Kip/Kiq/Kir);


%% -------------------------------
% 7. Succinate dehydrogenase
% SUC + COQ = QH2 + FUM

% load concentrations of reactants and products
A = SUC_x;
B = COQ_x;
P = QH2_x;
Q = FUM_x;

% dG and Keq values
dGr_sdho = dGf1(iQH2) + dGf1(iFUM) - dGf1(iSUC) - dGf1(iCOQ);
Keq_sdho = exp(-dGr_sdho/RT);
% Keq_sdh = Keq_sdho;
Keq_sdh = Keq_sdho*(P_x(iFUM)*P_x(iQH2))/(P_x(iSUC)*P_x(iCOQ));

% Km and Ki values (Molar)
KmA = 467e-6;
KmB = 480e-6;
KmP = 2.45e-6;
KmQ = 1200e-6;
Kia = 120e-6;
Kiq = 1275e-6;
KiOAA = 1.5e-6;
% % from Gopher and Gutman
% KaSUC = 800e-6;
% KaFUM = 6400e-6;
% from Kohn et al.
KaSUC = 450e-6;
KaFUM = 375e-6;

% inhibition coefficients
ai = (1+OAA_x/KiOAA+SUC_x/KaSUC+FUM_x/KaFUM)/(1+SUC_x/KaSUC+FUM_x/KaFUM);

% Vm values
Vmf = Vmax(isdh);
Vmr = Vmf/Keq_sdh*(KmP*Kiq/Kia/KmB);

% total reaction flux
J_sdh = (Vmf*Vmr*A*B - Vmf*Vmr*(P*Q/Keq_sdh)) / (Vmr*Kia*KmB*ai+Vmr*KmB*A ...
           +Vmr*KmA*ai*B+Vmf*KmQ*ai/Keq_sdh*P+Vmf*KmP/Keq_sdh*Q ...
           +Vmr*A*B+Vmf*KmQ/Kia/Keq_sdh*A*P+Vmr*KmA/Kiq*B*Q ...
           +Vmf/Keq_sdh*P*Q);


%% -------------------------------
% 8. Fumarase
% FUM (+ H2O) = MAL

% load concentrations of reactants and products
A = FUM_x;
P = MAL_x;

% dG and Keq values
dGr_fumo = dGf1(iMAL) - dGf1(iFUM) - dGf1(iH2O);
Keq_fumo = exp(-dGr_fumo/RT);
% Keq_fum = Keq_fumo;
Keq_fum = Keq_fumo*P_x(iMAL)/P_x(iFUM);

% Km and Ki values (Molar)
KmA = 44.7e-6;
KmP = 197.7e-6;
% KmA = 2.34e-6;
% KmP = 8e-6;
KiCIT = 3500e-6;
KiATP = 40e-6;
KiADP = 400e-6;
KiGTP = 80e-6;
KiGDP = 330e-6;

% inhibition coefficients
fATP_x = ATP_x * (1+H_x/Kh(iATP))/P_x(iATP);
fADP_x = ADP_x * (1+H_x/Kh(iADP))/P_x(iADP);
fGTP_x = GTP_x * (1+H_x/Kh(iGTP))/P_x(iGTP);
fGDP_x = GDP_x * (1+H_x/Kh(iGDP))/P_x(iGDP);
ai = 1+CIT_x/KiCIT+fATP_x/KiATP+fADP_x/KiADP+fGTP_x/KiGTP+fGDP_x/KiGDP;

% Vm values
Vmf = Vmax(ifum);
Vmr = Vmf/Keq_fum*(KmP/KmA);

% total reaction flux
J_fum =  (Vmf*Vmr*A - Vmf*Vmr*(P/Keq_fum))/(KmA*Vmr*ai+Vmr*A+Vmf/Keq_fum*P);


%% -------------------------------
% 9. Malate dehydrogenase
% NAD + MAL = OAA + NADH (+ H^+)

% load concentrations of reactants and products
A = NAD_x;
B = MAL_x;
P = OAA_x;
Q = NADH_x;

% dG and Keq values
dGr_mdho = dGf1(iOAA) + dGf1(iNADH) - dGf1(iNAD) - dGf1(iMAL);
Keq_mdho = exp(-dGr_mdho/RT);
% Keq_mdh = Keq_mdho*1/H_x;
Keq_mdh = Keq_mdho*1/H_x*(P_x(iOAA)*P_x(iNADH))/(P_x(iMAL)*P_x(iNAD));

% Km and Ki values (Molar)
KmA = 90.55e-6;
KmB = 250e-6;
KmP = 6.128e-6;
KmQ = 2.58e-6;
Kia = 279e-6;
Kib = 360e-6;
Kip = 5.5e-6;
Kiq = 3.18e-6;
% % from Kohn et al.
% KiATP = 709.3e-6;
% KiADP = 383.2e-6;
% KiAMP = 793.0e-6;
% from Oza and Shore
KiATP = 183.2e-6;
KiADP = 394.4e-6;
KiAMP = 420.0e-6;

% inhibition coefficients
fATP_x = ATP_x * (1+H_x/Kh(iATP))/P_x(iATP);
fADP_x = ADP_x * (1+H_x/Kh(iADP))/P_x(iADP);
fAMP_x = AMP_x * (1+H_x/Kh(iAMP))/P_x(iAMP);
ai = 1+fATP_x/KiATP+fADP_x/KiADP+fAMP_x/KiAMP;

% Vm values
Vmf = Vmax(imdh);
Vmr = Vmf/Keq_mdh*(Kiq*KmP/Kia/KmB);

% total reaction flux
J_mdh = (Vmf*Vmr*A*B - Vmf*Vmr*(P*Q/Keq_mdh)) / (Vmr*Kia*KmB*ai+Vmr*KmB*A ...
            +Vmr*KmA*ai*B+Vmf*KmQ*ai/Keq_mdh*P+Vmf*KmP/Keq_mdh*Q ...
            +Vmr*A*B+Vmf*KmQ/Kia/Keq_mdh*A*P+Vmf/Keq_mdh*P*Q ...
            +Vmr*KmA/Kiq*B*Q+Vmr/Kip*A*B*P+Vmf/Kib/Keq_mdh*B*P*Q);
            

%% -------------------------------
% 10. Nucleoside diphosphokinase
% GTP + ADP = GDP + ATP

% load concentrations of reactants and products
A = GTP_x;
B = ADP_x;
P = GDP_x;
Q = ATP_x;

% dG and Keq values
dGr_ndko = dGf1(iGDP) + dGf1(iATP) - dGf1(iGTP) - dGf1(iADP);
Keq_ndko = exp(-dGr_ndko/RT);
% Keq_ndk = Keq_ndko;
Keq_ndk = Keq_ndko*(P_x(iGDP)*P_x(iATP))/(P_x(iGTP)*P_x(iADP));

% Km and Ki values (Molar)
KmA = 111e-6;
KmB = 100e-6;
KmP = 260e-6;
KmQ = 278e-6;
Kia = 170e-6;
Kib = 143.6e-6;
Kip = 146.6e-6;
Kiq = 156.5e-6;
KiAMP = 650e-6;

% inhibition coefficients
fAMP_x = AMP_x * (1+H_x/Kh(iAMP))/P_x(iAMP);
ai = 1 + fAMP_x/KiAMP;

% Vm values
Vmf = Vmax(indk);
Vmr = Vmf/Keq_ndk*(KmQ*Kip/Kia/KmB);

% forward reaction flux
if (A > MinCon) && (B > MinCon) 
    J_ndk_f = Vmf*Vmr*A*B /ai / (Vmr*KmB*A+Vmr*KmA*B ...
            +Vmf*KmQ/Keq_ndk*P+Vmf*KmP/Keq_ndk*Q+Vmr*A*B ...
            +Vmf*KmQ/Kia/Keq_ndk*A*P+Vmf/Keq_ndk*P*Q ...
            +Vmr*KmA/Kiq*B*Q);
else
    J_ndk_f = 0;
end

% total reaction flux
if (P > MinCon) && (Q > MinCon)
    J_ndk = J_ndk_f - Vmf*Vmr*(P*Q/Keq_ndk)/ai / (Vmr*KmB*A+Vmr*KmA*B ...
            +Vmf*KmQ/Keq_ndk*P+Vmf*KmP/Keq_ndk*Q+Vmr*A*B ...
            +Vmf*KmQ/Kia/Keq_ndk*A*P+Vmf/Keq_ndk*P*Q ...
            +Vmr*KmA/Kiq*B*Q);
else
    J_ndk = J_ndk_f;
end


%% -------------------------------
% 11. Glutamate oxaloacetate transaminase (aspartate transaminase)
% ASP + AKG = OAA + GLU

% load concentrations of reactants and products
A = ASP_x;
B = AKG_x;
P = OAA_x;
Q = GLU_x;

% dG and Keq values
dGr_goto = dGf1(iOAA) + dGf1(iGLU) - dGf1(iASP) - dGf1(iAKG);
Keq_goto = exp(-dGr_goto/RT);
% Keq_got = Keq_goto;
Keq_got = Keq_goto*(P_x(iOAA)*P_x(iGLU))/(P_x(iASP)*P_x(iAKG));

% Km and Ki values (Molar)
KmA = 3900e-6;
KmB = 430e-6;
KmP = 88e-6;
KmQ = 8900e-6;
% Kia = 8400e-6;
% Kib = 50e-6;
% Kip = 710e-6;
% Kiq = 3480e-6;
Kia = 3480e-6;
Kib = 710e-6;
Kip = 50e-6;
Kiq = 8400e-6;
KiAKG = 16.6e-3;
ai = 1 + AKG_x/KiAKG;

% Vm values
Vmf = Vmax(igot);
Vmr = Vmf/Keq_got*(KmQ*Kip/Kia/KmB);
% Vmr1 = Vmf/Keq_got*(KmP*Kiq/Kib/KmA);

% forward reaction flux
if (A > MinCon) && (B > MinCon)
    J_got_f = Vmf*Vmr*A*B / (Vmr*KmB*A+Vmr*KmA*ai*B ...
            +Vmf*KmQ/Keq_got*ai*P+Vmf*KmP/Keq_got*Q+Vmr*A*B ...
            +Vmf*KmQ/Kia/Keq_got*A*P+Vmf/Keq_got*P*Q ...
            +Vmr*KmA/Kiq*B*Q); 
else 
    J_got_f = 0;
end

% total reaction flux
if (P > MinCon) && (Q > MinCon)
    J_got = J_got_f - Vmf*Vmr*(P*Q/Keq_got) / (Vmr*KmB*A+Vmr*KmA*ai*B ...
            +Vmf*KmQ/Keq_got*ai*P+Vmf*KmP/Keq_got*Q+Vmr*A*B ...
            +Vmf*KmQ/Kia/Keq_got*A*P+Vmf/Keq_got*P*Q ...
            +Vmr*KmA/Kiq*B*Q);
else
    J_got = J_got_f;
end

        
%% -------------------------------
% 12. Anti- and co-transporter fluxes of substrates involved in the TCA cycle
% Fluxes are defined to be positive when the first reactant(s) move from IM
% into mito. matrix

% --------------------------------
% (1) PYR^{-}-H^{+} co-transporter
PYR_i1 = PYR_i*1/P_i(iPYR);
PYR_x1 = PYR_x*1/P_x(iPYR);
J_PYR_H = x_PYR_H * (PYR_i1*H_i - PYR_x1*H_x);

% --------------------------------
% (2) GLU^{-}-H^{+} co-transporter
GLU_i1 = GLU_i*1/P_i(iGLU);
GLU_x1 = GLU_x*1/P_x(iGLU);
J_GLU_H = x_GLU_H * (GLU_i1*H_i - GLU_x1*H_x);

% --------------------------------
% (3) CIT^{2-}/MAL^{2-} anti-transporter
CIT_i1 = CIT_i*(H_i/Kh(iCIT))/P_i(iCIT);
CIT_x1 = CIT_x*(H_x/Kh(iCIT))/P_x(iCIT);
MAL_i1 = MAL_i*1/P_i(iMAL);
MAL_x1 = MAL_x*1/P_x(iMAL);
J_CIT_MAL = x_CIT_MAL * (CIT_i1*MAL_x1 - CIT_x1*MAL_i1);

% --------------------------------
% (4) AKG^{2-}/MAL^{2-} anti-transporter
KmMALi = 1.4e-3;
KmMALx = 0.7e-3;
KmAKGi = 0.3e-3;
KmAKGx = 0.17e-3;
J_AKG_MAL = x_AKG_MAL/(KmMALx*KmAKGi)*(MAL_x*AKG_i-MAL_i*AKG_x)/ ...
         (2 + MAL_i/KmMALi + MAL_x/KmMALx + AKG_i/KmAKGi + AKG_x/KmAKGx + ...
          MAL_i*AKG_x/(KmMALi*KmAKGx) + MAL_x*AKG_i/(KmMALx*KmAKGi));


% --------------------------------
% (5) MAL^{2-}/PI^{2-} anti-transporter
MAL_i1 = MAL_i*1/P_i(iMAL);
MAL_x1 = MAL_x*1/P_x(iMAL);
Pi_i1 = Pi_i*1/P_i(iPi);
Pi_x1 = Pi_x*1/P_x(iPi);
J_MAL_PI = x_MAL_PI * (MAL_i1*Pi_x1 - MAL_x1*Pi_i1);

% --------------------------------
% (6) ASP^{-}/H.GLU{0} antitransporter (with one negative charge
% translocated from IM into mito matrix)
Kiaspi = 28e-6;
Kiaspx = 2.8e-3;
Kiglui = 180e-6;
Kiglux = 1.6e-3;
pKa_gaa = 6.5;
Kh_gaa = 10^(-pKa_gaa);
Keq_gaa = exp(-F*dPsi/RT)*P_x(iASP)*P_i(iGLU)/P_i(iASP)/P_x(iGLU); % Modified, 06/11/08
m = 1.8;
J_ASP_GLU = x_ASP_GLU/(Keq_gaa*Kiaspi*Kiglux*Kh_gaa)*(Keq_gaa*ASP_i*GLU_x*H_x - ASP_x*GLU_i*H_i)/ ...
              (2*m + m*ASP_i/Kiaspi + ASP_i*GLU_x*H_x/(Kiaspi*Kiglux*Kh_gaa) + m*ASP_x*H_i/(Kiaspx*Kh_gaa) + ...
               ASP_x*GLU_i*H_i/(Kiaspx*Kiglui*Kh_gaa) + m*ASP_x/Kiaspx + m*ASP_i*H_x/(Kiaspi*Kh_gaa) + ...
               m*H_x/Kh_gaa + m*GLU_i*H_i/(Kiglui*Kh_gaa) + m*H_i/Kh_gaa + m*GLU_x*H_x/(Kiglux*Kh_gaa));


% --------------------------------
% (7) SUC^{2-}/MAL^{2-} anti-transporter
MAL_i1 = MAL_i*1/P_i(iMAL);
MAL_x1 = MAL_x*1/P_x(iMAL);
SUC_i1 = SUC_i*1/P_i(iSUC);
SUC_x1 = SUC_x*1/P_x(iSUC);
J_SUC_MAL = x_SUC_MAL * (SUC_i1*MAL_x1 - SUC_x1*MAL_i1);



%% -------------------------------
% 13. Passive permeation between cytoplasm/buffer and IM
% Fluxes are defined to be positive when the reactant moves from
% cytoplasm/buffer into IM

J_PYRt = gamma * x_PYRt * (PYR_c - PYR_i);
J_CITt = gamma * x_CITt * (CIT_c - CIT_i);
J_MALt = gamma * x_MALt * (MAL_c - MAL_i);
J_AKGt = gamma * x_AKGt * (AKG_c - AKG_i);
J_SUCt = gamma * x_SUCt * (SUC_c - SUC_i);
J_GLUt = gamma * x_GLUt * (GLU_c - GLU_i);
J_ASPt = gamma * x_ASPt * (ASP_c - ASP_i);
J_FUMt = gamma * x_FUMt * (FUM_c - FUM_i);
J_ICITt = gamma * x_ICITt * (ICIT_c - ICIT_i);


%% II. Flux expresssions in the oxidative phosphorylation 

%% ------------------------------- 
% Substrate/ion transport
% Transport of ADP, ATP, and Pi across outer membrane:
J_ADP   = gamma*x_A*(ADP_c-ADP_i);
J_ATP   = gamma*x_A*(ATP_c-ATP_i);
J_AMP   = gamma*x_A*(AMP_c-AMP_i);
J_Pi2   = gamma*x_Pi2*(Pi_c-Pi_i);    


%% -------------------------------
% 1. Complex I
% NADH_x + Q_x + 5H+_x <-> NAD+_x + QH2_x + 4H+_i + 4dPsi

% compute free energy of reaction from free energe of formation
dGr_C1o = dGf1(iNAD) + dGf1(iQH2) - dGf1(iNADH) - dGf1(iCOQ);
% compute Keq from combined free energy of reactions (including potential
% change due to charge translocation)
Keq_C1 = exp(-(dGr_C1o + 4*F*dPsi)/RT);
% compute Kapp from Kapp = Keq*H_x^n/H_i^m*Product(Poly)  (n,m are
% stochi.coefficients)
Kapp_C1 = Keq_C1*H_x^5/H_i^4;
J_C1 = x_C1*( Kapp_C1*NADH_x*COQ_x - NAD_x*QH2_x );

%% -------------------------------
% 2. Complex III
% QH2_x + 2cytoC(ox)3+_i + 2H+_x <-> Q_x + 2cytoC(red)2+_i + 4H+_i + 2dPsi

dGr_C3o = dGf1(iCOQ) + 2*dGf1(iCred) - dGf1(iQH2) - 2*dGf1(iCox);
Keq_C3 = exp(-(dGr_C3o + 2*F*dPsi)/RT);
Kapp_C3 = Keq_C3*H_x^2/H_i^4;
QH2_x = max(MinCon,QH2_x); COQ_x = max(MinCon,COQ_x);
J_C3 = x_C3*((1+Pi_x/k_Pi1)/(1+Pi_x/k_Pi2))*...
          (Kapp_C3^0.5*Cox_i*sqrt(QH2_x) - Cred_i*sqrt(COQ_x) );
      
%% -------------------------------
% 3. Complex IV
% 2cytoC(red)2+_i + 0.5O2_x + 4H+_x <-> 2cytoC(ox)3+_x + H2O_x + 2H+_i +
% 2dPsi

dGr_C4o = 2*dGf1(iCox) + dGf1(iH2O) - 2*dGf1(iCred) - 0.5*dGf1(iO2);
% 2 charges from translocation of proton, and the other 2 from cytoC
Keq_C4 = exp(-(dGr_C4o + 4*F*dPsi)/RT);
Kapp_C4 = Keq_C4*H_x^4/H_i^2;
O2 = max(O2,MinCon);
J_C4 = x_C4*(O2/(O2+k_O2))*exp(F*dPsi/RT)*(Cred_i/Ctot)*( Kapp_C4^0.5*Cred_i*(O2^0.25) - Cox_i );


%% -------------------------------
% 4. F1Fo-ATPase
% ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x

dGr_F1o = dGf1(iATP) + dGf1(iH2O) - dGf1(iADP) - dGf1(iPi);
Keq_F1 = exp(-(dGr_F1o-n_A*F*dPsi)/RT);
Kapp_F1 = Keq_F1*H_i^n_A/H_x^(n_A-1)*P_x(iATP)/(P_x(iADP)*P_x(iPi));
J_F1 = x_F1*(Kapp_F1*ADP_x*Pi_x - ATP_x);


%% -------------------------------
% 5. ANT
% ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

ADP_i1 = ADP_i/P_i(iADP); % ADP^3-
ATP_i1 = ATP_i/P_i(iATP); % ATP^4-
ADP_x1 = ADP_x/P_x(iADP); % ADP^3-
ATP_x1 = ATP_x/P_x(iATP); % ATP^4-

theta = 0.60;
Psi_i = theta*dPsi;
Psi_x = (theta-1)*dPsi;
if (ADP_i > 1e-9) && (ATP_x > 1e-9)
  J_ANT  = x_ANT*( ADP_i1/(ADP_i1+ATP_i1*exp(-F*Psi_i/RT)) - ADP_x1/(ADP_x1+ATP_x1*exp(-F*Psi_x/RT)) )*(ADP_i1/(ADP_i1+k_mADP));
else
  J_ANT  = 0;
end

%% -------------------------------
% 6. H+-Pi2 cotransporter

H2PIi1 = Pi_i*(H_i/Kh(iPi))/P_i(iPi);
H2PIx1 = Pi_x*(H_x/Kh(iPi))/P_x(iPi);
J_Pi1 = x_Pi1/k_PiH*(H_i*H2PIi1 - H_x*H2PIx1)/(1+H2PIi1/k_PiH)/(1+H2PIx1/k_PiH);


%% -------------------------------
% 7. H+ leaking
if abs(dPsi) > 1e-9
  J_Hle = x_Hle*dPsi*(H_i*exp(F*dPsi/RT)-H_x)/(exp(F*dPsi/RT)-1);
else
  J_Hle = x_Hle*RT*( H_i - H_x )/F;
end


%% -------------------------------
% 8. K+/H+ anti-porter
J_KH   = x_KH*( K_i*H_x - K_x*H_i);


%% -------------------------------
% 9. Rate of cytosolic ATP consumption by Hexokinase reaction
% GLC + ATP = ADP + G6P
% GLC^0 + ATP^4- = ADP^3- + G6P^2- + H^+
% Keq_HK = 880; % http://web.mit.edu/esgbio/www/eb/chem/solvingequilibria.html
dGr_HKo = dGf1(iG6P) + dGf1(iADP) - dGf1(iGLC) - dGf1(iATP);
Keq_HK = exp(-(dGr_HKo)/RT);
Kapp_HK = Keq_HK/H_c*P_c(iG6P)*P_c(iADP)/(P_c(iGLC)*P_c(iATP));

KmA = 1e-3; % MgATP
Kia = 1e-3;
KmB = 47e-6; % GLC
Kib = 47e-6;
KmP = 47e-6; % G6P
Kip = 47e-6;
KmQ = 1e-3; % MgADP
Kiq = 1e-3;
KiG6P = 10e-6; % inhibition constant of G6P
A = ATP_c*(Mg_c/Km(iATP))/P_c(iATP);
B = GLC_c;
P = G6P_c;
Q = ADP_c*(Mg_c/Km(iADP))/P_c(iADP);
Kapp_HK_m = Kapp_HK*(Km(iATP)*P_c(iATP))/(Km(iADP)*P_c(iADP));
J_HK = x_HK/(Kib*KmA)*(A*B - P*Q/Kapp_HK_m) ...
       /(1 + A/Kia + B/Kib + A*B/Kib/KmA + P/Kip + Q/Kiq + P*Q/Kiq/KmP + P*B/KiG6P/Kib);


%% ---------------------------------------
% 10. Creatine kinase reaction
% ADP3- + PCr- + H+ = ATP4- + Cr0
% set Cr and PCr concentrations according to experiments
if (ExpType == 1) || (ExpType == 2) || (ExpType == 5)
    J_CKe = 0;
elseif (ExpType == 7) %Muoio

x_CK = 1e7;
K_CK = exp(50.78/RT);                 
%CRtot = 42.7e-3; 
PCr_c = x(iPCr_c);
Cr_c  = 5e-3; %CRtot - PCr_c; % 5 mM for CK clamp %%%%5e-3; Muoio%%%%
ATP_c1 = ATP_c * 1/P_c(iATP); % Mg2+ unbound species
ADP_c1 = ADP_c * 1/P_c(iADP); % Mg2+ unbound species
J_CKe  = x_CK * (K_CK*ADP_c1*PCr_c*H_c - ATP_c1*Cr_c );
%     J_CKe = 0;

else 
x_CK = 1e7;
K_CK = exp(50.78/RT);                 
CRtot = 42.7e-3; 
PCr_c = x(iPCr_c);
Cr_c  =CRtot - PCr_c; % 5 mM for CK clamp 
ATP_c1 = ATP_c * 1/P_c(iATP); % Mg2+ unbound species
ADP_c1 = ADP_c * 1/P_c(iADP); % Mg2+ unbound species
J_CKe  = x_CK * (K_CK*ADP_c1*PCr_c*H_c - ATP_c1*Cr_c );
%     J_CKe = 0;

end


%% ------------------------------------------
% 11. Adenylate kinase reaction
% 2ADP3- = ATP4- + AMP2-
dGr_AKo = dGf1(iATP) + dGf1(iAMP) - 2*dGf1(iADP);
Keq_AK = exp(-dGr_AKo/RT);
Kapp_AKi = Keq_AK*P_i(iATP)*P_i(iAMP)/(P_i(iADP)*P_i(iADP)); % Added , 06/11/08
Kapp_AKc = Keq_AK*P_c(iATP)*P_c(iAMP)/(P_c(iADP)*P_c(iADP));
if (ExpType == 1) || (ExpType == 2) || (ExpType == 5) || (ExpType == 7) 
    J_AKi = 0;
    J_AKe = 0;
else
    x_AK = 1e7;
    J_AKi  = x_AK*( Kapp_AKi*ADP_i*ADP_i - AMP_i*ATP_i ); % Modified, 06/11/08
    J_AKe  = x_AK*( Kapp_AKc*ADP_c*ADP_c - AMP_c*ATP_c ); % Modified, 06/11/08
end


%% Computing time derivatives of state variables

%  (i) Matrix species and dPsi
f(idPsi)   = (1.48e5)*( 4*J_C1 + 2*J_C3 + 4*J_C4 - n_A*J_F1 - J_ANT - J_Hle + J_ASP_GLU);
f(iATP_x)  = (+J_ndk + J_F1 - J_ANT)/W_x; 
f(iADP_x)  = (-J_ndk - J_F1 + J_ANT)/W_x;
f(iAMP_x)  = (0)/W_x; 
f(iGTP_x)  = (+J_scoas - J_ndk)/W_x; 
f(iGDP_x)  = (-J_scoas + J_ndk)/W_x; 
f(iPi_x)   = (-J_scoas - J_F1 + J_Pi1  - J_MAL_PI)/W_x;
f(iNADH_x) = (+J_pdh + J_isod + J_akgd + J_mdh - J_C1)/W_x; 
f(iQH2_x)  = (+J_sdh + J_C1 - J_C3)/W_x; 
f(iPYR_x)  = (-J_pdh + J_PYR_H)/W_x; 
f(iACCOA_x) = (-J_cits + J_pdh)/W_x; 
f(iCIT_x)  = (+J_cits - J_acon + J_CIT_MAL)/W_x; 
f(iICIT_x) = (+J_acon - J_isod)/W_x; 
f(iAKG_x)  = (+J_isod - J_akgd - J_got + J_AKG_MAL)/W_x; 
f(iSCOA_x) = (+J_akgd - J_scoas)/W_x; 
f(iCOASH_x)  = (-J_pdh - J_akgd + J_scoas + J_cits)/W_x; 
f(iSUC_x)  = (+J_scoas - J_sdh + J_SUC_MAL)/W_x;
f(iFUM_x)  = (+J_sdh - J_fum )/W_x;  
f(iMAL_x)  = (+J_fum - J_mdh + J_MAL_PI - J_AKG_MAL - J_CIT_MAL - J_SUC_MAL)/W_x; 
f(iOAA_x)  = (-J_cits + J_mdh + J_got)/W_x;
f(iGLU_x)  = (+J_got + J_GLU_H - J_ASP_GLU)/W_x; 
f(iASP_x)  = (-J_got + J_ASP_GLU)/W_x; 
f(iO2_x)   = 0; 
f(iCO2tot_x) = 0;

%  (ii) IM space species
f(iCred_i) = (+2*J_C3 - 2*J_C4)/W_i; 
f(iATP_i) = (J_ATP + J_ANT + J_AKi)/W_i; 
f(iADP_i) = (J_ADP - J_ANT - 2*J_AKi)/W_i; 
f(iAMP_i) = (J_AMP + J_AKi)/W_i;   
f(iPi_i)  = (-J_Pi1 + J_Pi2 + J_MAL_PI)/W_i; 
f(iPYR_i) = (-J_PYR_H + J_PYRt)/W_i;
f(iCIT_i) = (-J_CIT_MAL + J_CITt)/W_i;
f(iICIT_i) = (J_ICITt)/W_i;
f(iAKG_i) = (-J_AKG_MAL + J_AKGt)/W_i;
f(iSUC_i) = (+ J_SUCt - J_SUC_MAL)/W_i;
f(iFUM_i) = (J_FUMt)/W_i;
f(iMAL_i) = (-J_MAL_PI + J_MALt + J_AKG_MAL + J_CIT_MAL + J_SUC_MAL)/W_i;
f(iGLU_i) = (-J_GLU_H + J_ASP_GLU + J_GLUt)/W_i;
f(iASP_i) = (-J_ASP_GLU + J_ASPt)/W_i;
f(iH_i)   = 0;
f(iMg_i)  = 0;
f(iK_i)   = 0;

%  (iii) Buffer species
% % d[ADP]/dt, d[ATP]/dt and d[PI]/dt in cytoplasm are assumed to be
% consumed by Hexokinase reaction
% W_c   [=] buffer water/mito || cyto water/cyto
% J_ATP [=] M/s/l mito 
% J_HK  [=] M/s/l buffer
% J_AtC, J_CKe [=] M/s/l cyto
% f [=] M/s/l cyto water
if (ExpType == 1) || (ExpType == 2) || (ExpType == 5) || (ExpType == 6)
    f(iPYR_c) = (-J_PYRt)/W_c;
    f(iCIT_c) = (-J_CITt)/W_c;
    f(iICIT_c) = (-J_ICITt)/W_c;
    f(iAKG_c) = (-J_AKGt)/W_c;
    f(iSUC_c) = (-J_SUCt)/W_c;
    f(iFUM_c) = (-J_FUMt)/W_c;
    f(iMAL_c) = (-J_MALt)/W_c;
    f(iGLU_c) = (-J_GLUt)/W_c;
    f(iASP_c) = (-J_ASPt)/W_c;
    f(iH_c)   = 0;
    f(iMg_c)  = 0;
    f(iK_c)   = 0;
    f(iPi_c)  = (-J_Pi2)/W_c; 
    f(iATP_c) = -J_HK/(W_c/(1+W_c))-J_ATP/W_c;  
    f(iADP_c) = +J_HK/(W_c/(1+W_c))-J_ADP/W_c; 
    f(iGLC_c) = -J_HK/(W_c/(1+W_c));  
    f(iG6P_c) = +J_HK/(W_c/(1+W_c));
    f(iPCr_c)  = 0;
    f(iAMP_c)  = 0;
elseif (ExpType == 7)
    f(iPYR_c) = (-J_PYRt)/W_c;
    f(iCIT_c) = (-J_CITt)/W_c;
    f(iICIT_c) = (-J_ICITt)/W_c;
    f(iAKG_c) = (-J_AKGt)/W_c;
    f(iSUC_c) = (-J_SUCt)/W_c;
    f(iFUM_c) = (-J_FUMt)/W_c;
    f(iMAL_c) = (-J_MALt)/W_c;
    f(iGLU_c) = (-J_GLUt)/W_c;
    f(iASP_c) = (-J_ASPt)/W_c;
    f(iH_c)   = 0;
    f(iMg_c)  = 0;
    f(iK_c)   = 0;
    f(iPi_c)  = (-J_Pi2+J_AtC)/W_c; %
    f(iADP_c) = (-J_ADP+J_AtC-J_CKe-2*J_AKe)/W_c;  %Changed 5/24 to account for Muoio experimental conditions LFM
    f(iGLC_c) =0; %-J_HK/W_c; %-J_HK/(W_c/(1+W_c));  %
    f(iG6P_c) = 0;%+J_HK/W_c; %0;%+J_HK/(W_c/(1+W_c)); %
    f(iATP_c) = 0;%(+J_ATP-J_AtC+J_CKe+J_AKe)/W_c;  Clamped
    f(iPCr_c) = 0; % Change
    f(iAMP_c) = (J_AMP+J_AKe)/W_c;
else
    f(iPYR_c) = (-Rm_cyto*J_PYRt)/W_c;
    f(iCIT_c) = (-Rm_cyto*J_CITt)/W_c;
    f(iICIT_c) = (-Rm_cyto*J_ICITt)/W_c;
    f(iAKG_c) = (-Rm_cyto*J_AKGt)/W_c;
    f(iSUC_c) = (-Rm_cyto*J_SUCt)/W_c;
    f(iFUM_c) = (-Rm_cyto*J_FUMt)/W_c;
    f(iMAL_c) = (-Rm_cyto*J_MALt)/W_c;
    f(iGLU_c) = (-Rm_cyto*J_GLUt)/W_c;
    f(iASP_c) = (-Rm_cyto*J_ASPt)/W_c;
    f(iH_c)   = 0;
    f(iMg_c)  = 0;
    f(iK_c)   = 0;
    f(iPi_c)  = (-Rm_cyto*J_Pi2+J_AtC)/W_c; 
    f(iATP_c) = (-Rm_cyto*J_ATP-J_AtC+J_CKe+J_AKe)/W_c;  
    f(iADP_c) = (-Rm_cyto*J_ADP+J_AtC-J_CKe-2*J_AKe)/W_c; 
    %f(iGLC_c) = 0;  
    %f(iG6P_c) = 0;  
    f(iPCr_c) = (-J_CKe)/W_c;
    f(iAMP_c) = (-Rm_cyto*J_AMP+J_AKe)/W_c;
    f(iAMP_c) = (-Rm_cyto*J_AMP+J_AKe)/W_c;
end


% Clamping the concentrations to generate st. st.
if (StateType == 0)
    f(iPYR_c) = 0;
    f(iMAL_c) = 0;
    f(iCIT_c) = 0;
    f(iICIT_c) = 0;
    f(iAKG_c) = 0;
    f(iSUC_c) = 0;
    f(iFUM_c) = 0;
end      

% Fixing [PYR] in the cytoplasm for in vivo experiments
if (ExpType == 3) || (ExpType == 4)|| (ExpType == 8)
    f(iPYR_c) = 0;
elseif (ExpType == 5) % for uncoupling
    f(iPYR_c) = 0;
    f(iMAL_c) = 0;
    f(iCIT_c) = 0;
    f(iICIT_c) = 0;
    f(iAKG_c) = 0;
    f(iSUC_c) = 0;
    f(iFUM_c) = 0;
end


%% Computing dH/dt, dMg/dt, and dK/dt
% Calculate dH/dt, dMg/dt, dK/dt according to Dan's strategy (available
% in the chapter "Biochemical Reaction Networks" in Beard and Qian's book)

% assume [H+], [Mg2+], [K+] constant in IM and cytoplasm/buffer space
% all the below caluclation if for dH_x/dt, dMg_x/dt, dK_x/dt

% take concentration of reactants located in the mito. matrix 
% from state varibles

xx(1:N_reactant) = 0;
xx(iH)      =   x(iH_x);
xx(iATP)    =   x(iATP_x);
xx(iADP)    =   x(iADP_x);
xx(iAMP)    =   x(iAMP_x);
xx(iGTP)    =   x(iGTP_x);
xx(iGDP)    =   x(iGDP_x);
xx(iPi)     =   x(iPi_x);
xx(iNADH)   =   x(iNADH_x);
xx(iNAD)    =   NADtot - x(iNADH_x);
xx(iQH2)    =   x(iQH2_x);
xx(iCOQ)    =   Qtot - x(iQH2_x);
xx(iPYR)    =   x(iPYR_x);
xx(iOAA)    =   x(iOAA_x);
xx(iACCOA)  =   x(iACCOA_x);
xx(iCIT)    =   x(iCIT_x);
xx(iICIT)   =   x(iICIT_x);
xx(iAKG)    =   x(iAKG_x);
xx(iSCOA)   =   x(iSCOA_x);
xx(iCOASH)  =   x(iCOASH_x);
xx(iSUC)    =   x(iSUC_x);
xx(iFUM)    =   x(iFUM_x);
xx(iMAL)    =   x(iMAL_x);
xx(iGLU)    =   x(iGLU_x);
xx(iASP)    =   x(iASP_x);
xx(iK)      =   x(iK_x);
xx(iMg)     =   x(iMg_x);
xx(iO2)     =   x(iO2_x);
xx(iFADH2)  =   FADtot/2; 
xx(iFAD)    =   FADtot - xx(iFADH2); 
xx(iCO2tot) =   x(iCO2tot_x); 

dxxdt(1:N_reactant) = 0;
% dxxdt(iH)   = f(iH_x);
dxxdt(iATP) = f(iATP_x);
dxxdt(iADP) = f(iADP_x);
dxxdt(iAMP) = f(iAMP_x);
dxxdt(iGTP) = f(iGTP_x);
dxxdt(iGDP) = f(iGDP_x);
dxxdt(iPi)  = f(iPi_x);
dxxdt(iNADH)= f(iNADH_x);
dxxdt(iNAD) = - f(iNADH_x);
dxxdt(iQH2) = f(iQH2_x);
dxxdt(iCOQ) = - f(iQH2_x);
dxxdt(iPYR) = f(iPYR_x);
dxxdt(iOAA) = f(iOAA_x);
dxxdt(iACCOA) = f(iACCOA_x);
dxxdt(iCIT) = f(iCIT_x);
dxxdt(iICIT)= f(iICIT_x);
dxxdt(iAKG) = f(iAKG_x);
dxxdt(iSCOA)= f(iSCOA_x);
dxxdt(iCOASH)= f(iCOASH_x);
dxxdt(iSUC) = f(iSUC_x);
dxxdt(iFUM) = f(iFUM_x);
dxxdt(iMAL) = f(iMAL_x);
dxxdt(iGLU) = f(iGLU_x);
dxxdt(iASP) = f(iASP_x);
% dxxdt(iK)   = f(iK_x);
% dxxdt(iMg)  = f(iMg_x);
dxxdt(iCred)= f(iCred_i);
dxxdt(iCox) = - f(iCred_i);
dxxdt(iO2)  = f(iO2_x);
dxxdt(iH2O) = 0;
dxxdt(iFADH2) = 0; 
dxxdt(iFAD) = - dxxdt(iFADH2); 
dxxdt(iCO2tot) = f(iCO2tot_x); 

% Partial Derivatives:
pHBpM = - sum( (H_x*xx./Kh)./(Km.*P_x.^2) );
pHBpK = - sum( (H_x*xx./Kh)./(Kk.*P_x.^2) );
pHBpH = + sum( (1+Mg_x./Km+K_x./Kk).*xx./(Kh.*P_x.^2) );
pMBpH = - sum( (Mg_x*xx./Km)./(Kh.*P_x.^2) );
pMBpK = - sum( (Mg_x*xx./Km)./(Kk.*P_x.^2) );
pMBpM = + sum( (1+H_x./Kh+K_x./Kk).*xx./(Km.*P_x.^2) );
pKBpH = - sum( (K_x*xx./Kk)./(Kh.*P_x.^2) );
pKBpM = - sum( (K_x*xx./Kk)./(Km.*P_x.^2) );
pKBpK = + sum( (1+H_x./Kh+Mg_x./Km).*xx./(Kk.*P_x.^2) );

% PHI's:
Phi_H = - sum( H_x*dxxdt./(Kh.*P_x) ) ...
        + 1*(-1*J_pdh + 2*J_cits + (-1)*J_akgd + J_scoas + J_mdh ... 
           + 1*(J_PYR_H + J_GLU_H + J_CIT_MAL  - J_ASP_GLU) ...
           -(4+1)*J_C1 - (4-2)*J_C3 - (2+2)*J_C4 + (n_A-1)*J_F1 ...
           + 2*J_Pi1 + J_Hle - 1*J_KH)/W_x;
Phi_M = - sum( Mg_x*dxxdt./(Km.*P_x) ) ;
Phi_K = - sum( K_x*dxxdt./(Kk.*P_x) ) + 1*J_KH/W_x ;


% alpha's:
aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;

% add additional buffer for [H+]
BX = 0.02; % M
K_BX = 1e-7; % M
aH = 1 + pHBpH + BX/K_BX/(1+H_x/K_BX)^2;

% Denominator:
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;

% Derivatives for H,Mg,K:
dH_xdt = ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
        
dMg_xdt = ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
        
dK_xdt = ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;      

f(iH_x) = dH_xdt;
f(iMg_x) = dMg_xdt;
f(iK_x) = dK_xdt;




%% Output the time derivative matrix
f = f';


%% -------------------------------
% end of function
% --------------------------------
