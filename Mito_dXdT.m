function [f] = Mito_dXdT(t,x,param,boxVmax,delta_sf)
% This function is used to calculate time derivatives of state
% variables of the cell-level cardaic energetics.

% Copyright: 
% Fan Wu and Daniel. A Beard, 2008
% Computational Engineering Group
% Biotechnology and Bioengineering Center
% Medical College of Wisconsin

% iflag: 0 - normal; 1 - no oxygen consumption; 
%       8 - fixed PI_c; 9 - fixed ADP_c;
%       10 - no CK activity; 11 -
%       no CK and fixed PI_c; 12 - no CK and fixed ADP_c;
%       100 - constant J_AtC = x_AtC;
%       101 - constant J_AtC = x_AtC and fixed [CrP]c and [Cr]_c
A=load('s_F');    
s_f=A.s_F;
% s_f(5)=delta_sf(1);
% s_f(6)=delta_sf(2);
% s_f(11)=delta_sf(3); 
% s_f(18)=delta_sf(4);
% s_f(19)=delta_sf(5);

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
% param(23) used by Kir1 in J_akgd
% param(24) not used 

% param(26) not used
% param(27) not used
% param(28) not used
% param(29) not used
% param(30) not used

% Parameters for oxidative phosphorylation
x_C1 = param(31);
x_C3 = param(32);
x_C4 = param(33);
x_F1 = param(34);
x_ANT = param(35);
x_PI1 = param(36);
k_PIH = param(37);
x_KH = param(38);
x_Hle = param(39);
k_PI1 = param(40);
k_PI2 = param(41);
ref_PYR = param(42); 


% Permeability coefficients for TCA intermediates
p_TI     = 85*1; % assumed to be equal to x_A for nucleotides (micro sec^-1)
x_PYRt   = p_TI;
x_GLUt   = p_TI;
x_ASPt   = p_TI;
x_CITt   = p_TI;
x_ICITt  = p_TI;
x_AKGt   = p_TI;
x_FUMt   = 0*p_TI;
x_SUCt   = p_TI;
x_MALt   = p_TI;

%% define indicis for all state variables
% (i) oxygen
iPO2         = 1;
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
iCr_c       = 62;
%%  Beta Oxidation indices
iC16Carn_cy     = 63;   % C16 AcylCarnitine cytosol
iC16Carn_m      = 64;   % C16 AcylCarnitine Matrix
iC16CoA_m       = 65;   % C16 AcylCoA Matrix
iC16EnoylCoA_m  = 66;   % C16 EnoylCoA Matrix
iC16OHCoA_m     = 67;   % C16 HydroxoxyacylCoA Matrix
iC16KetoCoA_m   = 68;   % C16 KetoacylCoA Matrix

iC14Carn_cy     = 69;   % C14 AcylCarnitine cytosol
iC14Carn_m      = 70;   % C14 AcylCarnitine Matrix
iC14CoA_m       = 71;   % C14 AcylCoA Matrix
iC14EnoylCoA_m  = 72;   % C14 EnoylCoA Matrix
iC14OHCoA_m     = 73;   % C14 HydroxoxyacylCoA Matrix
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
iAcetylCoAMAT   = 105; 
iFADH_m         = 106;
iNADHm          = 107; 
iCoAMAT         = 108; 
iC16AcylCoACYT  = 109;
iCARN_c         = 110;
iCARN_i         = 111;


C16Carn_cy      = x(63);   % C16 AcylCarnitine cytosol
C16Carn_m       = x(64);   % C16 AcylCarnitine Matrix
C16CoA_m        = x(65);   % C16 AcylCoA Matrix
C16EnoylCoA_m   = x(66);   % C16 EnoylCoA Matrix
C16OHCoA_m      = x(67);   % C16 HydroxoxyacylCoA Matrix
C16KetoCoA_m    = x(68);   % C16 KetoacylCoA Matrix

C14Carn_cy      = x(69);   % C14 AcylCarnitine cytosol
C14Carn_m       = x(70);   % C14 AcylCarnitine Matrix
C14CoA_m        = x(71);   % C14 AcylCoA Matrix
C14EnoylCoA_m   = x(72);   % C14 EnoylCoA Matrix
C14OHCoA_m      = x(73);   % C14 HydroxoxyacylCoA Matrix
C14KetoCoA_m    = x(74);

C12Carn_cy      = x(75); 
C12Carn_m       = x(76);
C12CoA_m        = x(77);
C12EnoylCoA_m   = x(78); 
C12OHCoA_m      = x(79);    
C12KetoCoA_m    = x(80); 

C10Carn_cy      = x(81); 
C10Carn_m       = x(82);   
C10CoA_m        = x(83);
C10EnoylCoA_m   = x(84);
C10OHCoA_m      = x(85);  
C10KetoCoA_m    = x(86); 

C8Carn_cy       = x(87); 
C8Carn_m        = x(88); 
C8CoA_m         = x(89); 
C8EnoylCoA_m    = x(90);
C8OHCoA_m       = x(91);
C8KetoCoA_m     = x(92); 

C6Carn_cy       = x(93); 
C6Carn_m        = x(94);     
C6CoA_m         = x(95);
C6EnoylCoA_m    = x(96);
C6OHCoA_m       = x(97);
C6KetoCoA_m     = x(98); 

C4Carn_cy       = x(99);  
C4Carn_m        = x(100);    
C4CoA_m         = x(101);
C4EnoylCoA_m    = x(102);
C4OHCoA_m       = x(103);
C4KetoCoA_m     = x(104); 

FADH_m          = x(106); 
  
C16AcylCoACYT   = x(109);

const_species_CarCYT = x(110); 
const_species_CarCYTi = x(111);


%% Defining indices of reactants
iH      =   1;
iATP    =   2;
iADP    =   3;
iAMP    =   4;
iGTP    =   5;
iGDP    =   6;
iPI     =   7;
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
RT = 8.314*(37+273.15)/1e3;          % kJ  mol^{-1}
F = 0.096484;                   % kJ mol^{-1} mV^{-1}

%  (ii) Subcelular volumes and water spaces
Vmito   = 0.001;                 % (ml mito / ml cuvette) %DAB 2/11/202
Vbuffer = 1-Vmito;               % (ml buffer / ml cuvette) %DAB 2/11/2020
Rm_buffer = Vmito / Vbuffer;     % Volume ratio mito volume / buffer volume

W_c = 0.807*1.044;              % cytosol water space (ml water per ml cytosol) [VB2002]
W_m = 0.664*1.09;               % mitochondrial water space (ml water per ml mito) [VB2002]
W_x = 0.9*W_m;                  % Matrix water space  (90% is matrix)
W_i = 0.1*W_m;                  % IM water space	(10% is intermembranee-space)

%  (iii) Pooled concentrations
SF=10^6; % uM -> Mf
Ctot   = 2.70e-3;               % M; total cytoC
Qtot   = 1.35e-3;               % M; total Q + QH2
NADtot = 250/SF; 
FADtot = 0.1e-3;                % M; total FADred+FADox (from Feng Yang)

%  (iv) Ox Phos Model Parameters
n_A     = 3.0;  % unitless
k_O2   = 1.2e-4;
% k_mADP = 3.5e-6;

%  (v) Outer membrane transport parameters
x_A    = 85;        % micron sec^{-1}
x_PI2  = 327;       % micron sec^{-1}
gamma  = 5.99;      % mito membrane area per cell volume micron^{-1}

% % (vi) Oxygen solubility and related parameters
a_3 = 1.74e-6;                    % oxygen solubility in cell
% CMb = 200e-6;                   % Myoglobin concentration (moles per liter cell volume)
% P50 = 2.39;                     % oxy-myoglobin 1/2 saturation

%% Loading values of state variables
MinCon = 1e-32;

% (i) Matrix species and dPsi
PO2         = x(iPO2);
dPsi        = x(idPsi);
ATP_x       = x(iATP_x);
ADP_x       = x(iADP_x);
AMP_x       = x(iAMP_x);
GTP_x       = x(iGTP_x);
GDP_x       = x(iGDP_x);
PI_x        = x(iPI_x);
NADH_x      = x(iNADH_x);
%
NADHm = NADH_x;    %%  

QH2_x       = x(iQH2_x);
PYR_x       = x(iPYR_x);
OAA_x       = x(iOAA_x);
ACCOA_x     = x(iACCOA_x);  
%
AcetylCoAMAT = ACCOA_x;

CIT_x       = x(iCIT_x);
ICIT_x      = x(iICIT_x);
AKG_x       = x(iAKG_x);
SCOA_x      = x(iSCOA_x);

const_species_CoAMATt=5000.0/SF; %M
CoAMAT=const_species_CoAMATt-(C16CoA_m+C16EnoylCoA_m+C16OHCoA_m+C16KetoCoA_m+C14CoA_m+C14EnoylCoA_m+C14OHCoA_m+C14KetoCoA_m+C12CoA_m+C12EnoylCoA_m+C12OHCoA_m+C12KetoCoA_m+C10CoA_m+C10EnoylCoA_m+C10OHCoA_m+C10KetoCoA_m+C8CoA_m+C8EnoylCoA_m+C8OHCoA_m+C8KetoCoA_m+C6CoA_m+C6EnoylCoA_m+C6OHCoA_m+C6KetoCoA_m+C4CoA_m+C4EnoylCoA_m+C4OHCoA_m+C4KetoCoA_m+AcetylCoAMAT+SCOA_x);

COASH_x     = CoAMAT; 

SUC_x       = x(iSUC_x);
FUM_x       = x(iFUM_x);
MAL_x       = x(iMAL_x);
GLU_x       = x(iGLU_x);
ASP_x       = x(iASP_x);
H_x         = x(iH_x); 
K_x         = x(iK_x);
Mg_x        = x(iMg_x);
CO2tot      = x(iCO2tot_x);

% (ii) IM space species
Cred_i      = max(0,x(iCred_i));
ATP_i       = x(iATP_i);
ADP_i       = x(iADP_i);
AMP_i       = x(iAMP_i);
PI_i        = x(iPI_i);
PYR_i       = x(iPYR_i);
CIT_i       = x(iCIT_i);
AKG_i       = x(iAKG_i);
SUC_i       = x(iSUC_i);
MAL_i       = x(iMAL_i);
GLU_i       = x(iGLU_i);
ASP_i       = x(iASP_i);
% H_i         = x(iH_i); 
% Mg_i        = x(iMg_i);
% K_i         = x(iK_i);
FUM_i       = x(iFUM_i);
ICIT_i      = x(iICIT_i);

% (iii) Cytoplasmic species
ATP_c       = x(iATP_c);
ADP_c       = x(iADP_c);
PI_c        = x(iPI_c);
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
AMP_c       = x(iAMP_c);

% (iv) Other concentrations computed from the state variables:
NAD_x  = NADtot - NADH_x;
COQ_x  = Qtot - QH2_x;
Cox_i  = Ctot - Cred_i;

% (v) set the H+, Mg2+, and K+ are permeable for outer mito membrane
H_i = H_c;
Mg_i = Mg_c;
K_i = K_c;

% Oxygen concentrations:
CfO2 = a_3*PO2;

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
dGf1(iPI) = -1098.27;        % HPO42-
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
Kh(iPI) = 10^(-6.71); Km(iPI) = 10^(-1.69); Kk(iPI) = 10^(+0.0074);
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
P_x(iPI) = 1 + H_x/Kh(iPI) + Mg_x/Km(iPI) + K_x/Kk(iPI); % add K-bound item, 06/10/08
P_c(iPI) = 1 + H_c/Kh(iPI) + Mg_c/Km(iPI) + K_c/Kk(iPI); % add K-bound item, 06/10/08
P_i(iPI) = 1 + H_i/Kh(iPI) + Mg_i/Km(iPI) + K_i/Kk(iPI); % add K-bound item, 06/10/08
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

% PDH flux changed by DAB 3/5/2020

%PDH_matrix
a = PYR_x;
b = COASH_x;
c = NAD_x;
p = CO2tot;
q = ACCOA_x;

% dG and Keq vlaues
dGr_pdho = dGf1(iCO2tot) + dGf1(iACCOA) + dGf1(iNADH) ...
          - dGf1(iPYR) - dGf1(iCOASH) - dGf1(iNAD) - dGf1(iH2O);
Keq_pdho = exp(-dGr_pdho/RT);
Keq_pdh = Keq_pdho*(1/H_x)*(P_x(iCO2tot)*P_x(iACCOA)*P_x(iNADH)) ...
                           /(P_x(iPYR)*P_x(iCOASH)*P_x(iNAD));
                       
% get free NADH
Kn_NADH=0.3e-3; % NADH binding dissociation coefficient
Xcp0_NADH=3.5e-3; % NADH binding capacity
r=(NADH_x-Kn_NADH-Xcp0_NADH+sqrt((NADH_x-Kn_NADH-Xcp0_NADH).^2+4*Kn_NADH.*NADH_x))/2;
KmA=38.3e-6;
KmB=9.9e-6;
KmC=60.7e-6;
KiACCOA=40.2e-6;
KiNADH=40.0e-6;

% PDH Activation (DAB):
rat=ACCOA_x/COASH_x*ATP_x/ADP_x*NADH_x/NAD_x;
if(rat>0.9999)&&(rat<1.0001)
  aa=0.5;
else
  aa=(3*rat-1-sqrt(9*rat^2-14*rat+9))/(4*rat-4);
end
Vmf = aa*Vmax(ipdh);  

%Inhibitionconstants
ai1=1+ACCOA_x/KiACCOA;
ai2=1+r/KiNADH;

J_pdh=Vmf*(a*b*c-p*q*r/Keq_pdh)/(KmC*ai2*a*b+KmB*ai1*a*c+KmA*b*c+a*b*c);

%[J_pdh]=PDHphosdephos(ACCOA_x,ATP_x,ADP_x,NADH_x,CoAMAT,NAD_x,PYR_x,CO2tot,Vmf,Keq_pdh,ai1,ai2,ref_PYR,Varstruc);

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
% Keq_cits = Keq_citso*(1/H_x^2);
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
C = PI_x;
P = COASH_x;
Q = SUC_x;
R = GTP_x;

    
% dG and Keq values
dGr_scoaso = dGf1(iCOASH) + dGf1(iSUC) + dGf1(iGTP) ...
            - dGf1(iGDP) - dGf1(iSCOA) - dGf1(iPI);
Keq_scoaso = exp(-dGr_scoaso/RT);
% Keq_scoas = Keq_scoaso*(1/H_x);
Keq_scoas = Keq_scoaso*(1/H_x)*(P_x(iCOASH)*P_x(iSUC)*P_x(iGTP)) ...
                           /(P_x(iSCOA)*P_x(iPI)*P_x(iGDP));

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
Keq_got = Keq_goto*(P_x(iOAA)*P_x(iGLU))/((1 + H_x/Kh(iASP))*P_x(iAKG));

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
PI_i1 = PI_i*1/P_i(iPI);
PI_x1 = PI_x*1/P_x(iPI);
J_MAL_PI = x_MAL_PI * (MAL_i1*PI_x1 - MAL_x1*PI_i1);

% --------------------------------
% (6) ASP^{-}/H.GLU{0} antitransporter (with one negative charge
% translocated from IM into mito matrix)
Kiaspi = 28e-6;
Kiaspx = 2.8e-3;
Kiglui = 180e-6;
Kiglux = 1.6e-3;
pKa_gaa = 6.5;
Kh_gaa = 10^(-pKa_gaa);
Keq_gaa = exp(-F*dPsi/RT)*P_x(iASP)*P_c(iGLU)/P_c(iASP)/P_x(iGLU);
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
% Transport of ADP, ATP, and PI across outer membrane:
J_ADP   = gamma*x_A*(ADP_c-ADP_i);
J_ATP   = gamma*x_A*(ATP_c-ATP_i);
J_AMP   = gamma*x_A*(AMP_c-AMP_i);
J_PI2   = gamma*x_PI2*(PI_c-PI_i);    
J_CARN  = gamma*p_TI*(const_species_CarCYT-const_species_CarCYTi); %% LFM 4/3/ 2020 add permeability of Carnitine

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
J_C3 = x_C3*((1+PI_x/k_PI1)/(1+PI_x/k_PI2))*...
          (Kapp_C3^0.5*Cox_i*sqrt(QH2_x) - Cred_i*sqrt(COQ_x) );

%% -------------------------------
% 3. Complex IV
% 2cytoC(red)2+_i + 0.5O2_x + 4H+_x <-> 2cytoC(ox)3+_x + H2O_x + 2H+_i +
% 2dPsi

dGr_C4o = 2*dGf1(iCox) + dGf1(iH2O) - 2*dGf1(iCred) - 0.5*dGf1(iO2);
% 2 charges from translocation of proton, and the other 2 from cytoC
Keq_C4 = exp(-(dGr_C4o + 4*F*dPsi)/RT);
Kapp_C4 = Keq_C4*H_x^4/H_i^2;
J_C4 = x_C4*(CfO2/(CfO2+k_O2))*exp(F*dPsi/RT)*(Cred_i/Ctot)*( Kapp_C4^0.5*Cred_i*(CfO2^0.25) - Cox_i );


%% -------------------------------
% 4. F1Fo-ATPase
% ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x

dGr_F1o = dGf1(iATP) + dGf1(iH2O) - dGf1(iADP) - dGf1(iPI);
Keq_F1 = exp(-(dGr_F1o-n_A*F*dPsi)/RT);
Kapp_F1 = Keq_F1*H_i^n_A/H_x^(n_A-1)*P_x(iATP)/(P_x(iADP)*P_x(iPI));
J_F1 = x_F1*(Kapp_F1*ADP_x*PI_x - ATP_x);


%% -------------------------------
% 5. ANT
% ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

ADP_i1 = ADP_i/P_i(iADP); % ADP^3-
ATP_i1 = ATP_i/P_i(iATP); % ATP^4-
ADP_x1 = ADP_x/P_x(iADP); % ADP^3-
ATP_x1 = ATP_x/P_x(iATP); % ATP^4-

del_D = 0.0167;
del_T = 0.0699;
k2_ANT = 9.54/60; % = 1.59e-1
k3_ANT = 30.05/60; % = 5.01e-1
K_D_o_ANT = 38.89e-6; 
K_T_o_ANT = 56.05e-6;
A = +0.2829;
B = -0.2086;
C = +0.2372;
fi = F*dPsi/RT;
k2_ANT_fi = k2_ANT*exp((A*(-3)+B*(-4)+C)*fi);
k3_ANT_fi = k3_ANT*exp((A*(-4)+B*(-3)+C)*fi);

K_D_o_ANT_fi = K_D_o_ANT*exp(3*del_D*fi);
K_T_o_ANT_fi = K_T_o_ANT*exp(4*del_T*fi);
q = k3_ANT_fi*K_D_o_ANT_fi*exp(fi)/(k2_ANT_fi*K_T_o_ANT_fi);
term2 = k2_ANT_fi*ATP_x1.*ADP_i1*q/K_D_o_ANT_fi ;
term3 = k3_ANT_fi.*ADP_x1.*ATP_i1/K_T_o_ANT_fi;
num = term2 - term3;
den = (1 + ATP_i1/K_T_o_ANT_fi + ADP_i1/K_D_o_ANT_fi)*(ADP_x1 + ATP_x1*q);
J_ANT = x_ANT/7.2679e-003*(0.70e-1)*num/den; % x_ANT'(in the paper) = x_ANT/7.2679e-003*(0.70e-1)

%% -------------------------------
% 6. H+-PI2 cotransporter

H2PIi1 = PI_i*(H_i/Kh(iPI))/P_i(iPI);
H2PIx1 = PI_x*(H_x/Kh(iPI))/P_x(iPI);
J_PI1 = x_PI1/k_PIH*(H_i*H2PIi1 - H_x*H2PIx1)/(1+H2PIi1/k_PIH)/(1+H2PIx1/k_PIH);


%% -------------------------------
% 7. H+ leak
if abs(dPsi) > 1e-9
  J_Hle = x_Hle*dPsi*(H_i*exp(F*dPsi/RT)-H_x)/(exp(F*dPsi/RT)-1);
else
  J_Hle = x_Hle*RT*( H_i - H_x )/F;
end


%% -------------------------------
% 8. K+/H+ anti-porter
J_KH   = x_KH*( K_i*H_x - K_x*H_i);


% 
% %% ---------------------------------------
% % 10. Creatine kinase reaction
% % ADP3- + PCr2- + H+ = ATP4- + Cr0
% % set Cr and PCr concentrations according to experiments
% 
% x_CK = 1e7;
% % K_CK = exp(50.78/RT); 
% K_CK = 7.408e8; 
% Kapp_CK = K_CK*H_c*P_c(iATP)*P_c(iCr)/P_c(iADP)/P_c(iPCr);
% PCr_c = x(iPCr_c);
% Cr_c   = x(iCr_c); %CRtot - PCr_c;
% J_CKe  = x_CK * (Kapp_CK*ADP_c*PCr_c - ATP_c*Cr_c );
% 
% 
% %% ------------------------------------------
% % 11. Adenylate kinase reaction
% % 2ADP3- = ATP4- + AMP2-
% % dGr_AKo = dGf1(iATP) + dGf1(iAMP) - 2*dGf1(iADP);
% % Keq_AK = exp(-dGr_AKo/RT);
% Keq_AK = 3.97e-1;
% Kapp_AK = Keq_AK*P_c(iATP)*P_c(iAMP)/(P_c(iADP)*P_c(iADP));
% x_AK = 1e7;
% J_AKi  = 0*x_AK*( Kapp_AK*ADP_i*ADP_i - AMP_i*ATP_i );
% J_AKe  = x_AK*( Kapp_AK*ADP_c*ADP_c - AMP_c*ATP_c );
% 
% %% -------------------------------------------
% % 12. ATPase flux
% % ATP4- + H2O = ADP3- + PI2- + H+
% % x_AtC [=] mmol s^{-1} (l cell}^{-1}
% % J_AtC [=] mmol s^{-1} (l cyto}^{-1}
% % dGr_AtCo = dGf1(iADP) + dGf1(iPI) - dGf1(iATP) - dGf1(iH2O) ;
% % Keq_AtC = exp(-dGr_AtCo/RT);
% Keq_AtC = 1.16e-1; % updated based on the adjusted thermo data at 310.15K
% Kapp_AtC = Keq_AtC*1/H_c*P_c(iADP)*P_c(iPI)/P_c(iATP);
% dGr_AtCapp = -RT*log(Kapp_AtC);
% dGr_AtC = dGr_AtCapp + RT*log(PI_c*ADP_c/ATP_c);
% J_AtC  = (1/Rc_cell)*x_AtC*(1) / ( 1 + (50*PI_c*ADP_c/ATP_c) ); 


%%  ---------------------------------------------------------------------------
%%  Beta Oxidation model (M/sec)
% This portion of the model is derived from doi:10.1371/journal.pcbi.1003186
% NADH, FADH2 and Acetyl-CoA were treated as fixed external parameters
% (conc. of each < K1xsink parameter)
% Beta-oxidation parameters are M/sec or umol/s/mg protein.
% SF_enzymeCn is the specificity factor that determines the enzyme activity
% for the substrate with a specific chain length as a percentage of Vmax
% 
% Changes from base model: CoAMAT = COASH, AcetylCoAMAT =  ACCOA_x, NADHm =
% NADH_x, J(iCoAMATO) = 0  
% Step 1. Add in flux 
% Step 2. Set NADH_m to NADH_x/1000
% Step 3. CoA MAT replaced with COASH_x  
% Step 4. AcetylCoA Mat to ACCOA_x/1000
% Step 5. Remove Sink from ACCOA and NADH flux equations 
% Step 6. Add fluxes to AccoA and NADH *1000 (boxnadh, boxaccoa)
%%  ---------------------------------------------------------------------------
%% Conversion paramters
rho_m = 3.6697e-6;  % (l mito) (mg protein)^{-1} 
protein = 0.5*2; % mg/mL * 2mL % 1 mg, 

SF=10^6; % uM -> M
SF2=(10^6);%/60; M/s
SF3 =60; % s/min
 
const_species_CoACYT=0; 	   %M 
const_species_MalCoACYT=0.0;

const_species_CarMAT=950.0/SF; %M

const_species_FADtMAT=0.77/SF; %M
const_species_NADtMAT=250.0/SF;


%% Functions

function z=pow(x,y)
z=x^y;
end

function z=root(x,y);z=y^(1/x);
end

function z = piecewise(varargin)
	numArgs = nargin;
	result = 0;
	foundResult = 0;
	for k=1:2: numArgs-1
		if varargin{k+1} == 1
			result = varargin{k};
			foundResult = 1;
			break;
        end
	end
	if foundResult == 0
		result = varargin{numArgs};
    end
z = result;
end

function [z]=function_CPT1(sf,V,Kms1,Kms2,Kmp1,Kmp2,Ki1,Keq,S1,S2,P1,P2,I1,n)
z=(sf*V*(S1*S2/(Kms1*Kms2)-P1*P2/(Kms1*Kms2*Keq))/((1+S1/Kms1+P1/Kmp1+(I1/Ki1)^n)*(1+S2/Kms2+P2/Kmp2)));
end

function z=function_CACT(Vf,Vr,Kms1,Kms2,Kmp1,Kmp2,Kis1,Kip2,Keq,S1,S2,P1,P2)
z=(Vf*(S1*S2-P1*P2/Keq)/(S1*S2+Kms2*S1+Kms1*S2*(1+P2/Kip2)+Vf/(Vr*Keq)*(Kmp2*P1*(1+S1/Kis1)+P2*(Kmp1+P1))));
end

function z=function_CPT2(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms6,Kms7,Kms8,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Kmp7,Kmp8,Keq,S1,S2,S3,S4,S5,S6,S7,S8,P1,P2,P3,P4,P5,P6,P7,P8);
z=(sf*V*(S1*S8/(Kms1*Kms8)-P1*P8/(Kms1*Kms8*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5+S6/Kms6+P6/Kmp6+S7/Kms7+P7/Kmp7)*(1+S8/Kms8+P8/Kmp8)));
end

function z=function_VLCAD(sf,V,Kms1,Kms2,Kms3,Kms4,Kmp1,Kmp2,Kmp3,Kmp4,Keq,S1,S2,S3,S4,P1,P2,P3,P4);
z=(sf*V*(S1*(S4-P4)/(Kms1*Kms4)-P1*P4/(Kms1*Kms4*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3)*(1+(S4-P4)/Kms4+P4/Kmp4)));
end

function z=function_LCAD(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms6,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Keq,S1,S2,S3,S4,S5,S6,P1,P2,P3,P4,P5,P6);
z=(sf*V*(S1*(S6-P6)/(Kms1*Kms6)-P1*P6/(Kms1*Kms6*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5)*(1+(S6-P6)/Kms6+P6/Kmp6)));
end

function z=function_MCAD(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms6,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Keq,S1,S2,S3,S4,S5,S6,P1,P2,P3,P4,P5,P6);
z=(sf*V*(S1*(S6-P6)/(Kms1*Kms6)-P1*P6/(Kms1*Kms6*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5)*(1+(S6-P6)/Kms6+P6/Kmp6)));
end

function z=function_SCAD(sf,V,Kms1,Kms2,Kms3,Kmp1,Kmp2,Kmp3,Keq,S1,S2,S3,P1,P2,P3);
z=(sf*V*(S1*(S3-P3)/(Kms1*Kms3)-P1*P3/(Kms1*Kms3*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2)*(1+(S3-P3)/Kms3+P3/Kmp3)));
end

function z=function_CROT(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms6,Kms7,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Kmp7,Ki1,Keq,S1,S2,S3,S4,S5,S6,S7,P1,P2,P3,P4,P5,P6,P7,I1);
z=(sf*V*(S1/Kms1-P1/(Kms1*Keq))/(1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5+S6/Kms6+P6/Kmp6+S7/Kms7+P7/Kmp7+I1/Ki1));
end

function z=function_MSCHAD(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms6,Kms7,Kms8,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Kmp7,Kmp8,Keq,S1,S2,S3,S4,S5,S6,S7,S8,P1,P2,P3,P4,P5,P6,P7,P8);
z=(sf*V*(S1*(S8-P8)/(Kms1*Kms8)-P1*P8/(Kms1*Kms8*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5+S6/Kms6+P6/Kmp6+S7/Kms7+P7/Kmp7)*(1+(S8-P8)/Kms8+P8/Kmp8)));
end

function z=function_MCKATA(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms6,Kms7,Kms8,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Kmp7,Kmp8,Keq,S1,S2,S3,S4,S5,S6,S7,S8,P1,P2,P3,P4,P5,P6,P7,P8);
 z=(sf*V*(S1*S8/(Kms1*Kms8)-P1*P8/(Kms1*Kms8*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5+S6/Kms6+P6/Kmp6+S7/Kms7+P7/Kmp7+P8/Kmp8)*(1+S8/Kms8+P8/Kmp8)));
end

function z=function_MCKATB(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms6,Kms7,Kms8,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Kmp7,Kmp8,Keq,S1,S2,S3,S4,S5,S6,S7,S8,P1,P2,P3,P4,P5,P6,P7,P8);
z=(sf*V*(S1*S8/(Kms1*Kms8)-P8*P8/(Kms1*Kms8*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5+S6/Kms6+P6/Kmp6+S7/Kms7+P7/Kmp7+P8/Kmp8)*(1+S8/Kms8+P8/Kmp8)));
end

function z=function_MTP(sf,V,Kms1,Kms2,Kms3,Kms4,Kms5,Kms7,Kms8,Kmp1,Kmp2,Kmp3,Kmp4,Kmp5,Kmp6,Kmp7,Kmp8,Ki1,Keq,S1,S2,S3,S4,S5,S7,S8,P1,P2,P3,P4,P5,P6,P7,P8,I1);
z=(sf*V*(S1*(S7-P7)*S8/(Kms1*Kms7*Kms8)-P1*P7*P8/(Kms1*Kms7*Kms8*Keq))/((1+S1/Kms1+P1/Kmp1+S2/Kms2+P2/Kmp2+S3/Kms3+P3/Kmp3+S4/Kms4+P4/Kmp4+S5/Kms5+P5/Kmp5+P6/Kmp6+I1/Ki1)*(1+(S7-P7)/Kms7+P7/Kmp7)*(1+S8/Kms8+P8/Kmp8)));
end

function z=function_RES(Ks,S,K1); 
z=(Ks*(S-K1));
end

%% Beta oxidation Paramters	
    compartment_VCYT = 10E-2; % L cyto/mg protein      
	compartment_VMAT = 1.8E-6;  
 
    VMAT = 1;
	global_par_Vfcact=boxVmax(1);       
	global_par_Vrcact=boxVmax(2);       
	global_par_KmcactCarMAT=130.0/SF;   % M
	global_par_KmcactCarCYT=130.0/SF;
	global_par_KicactCarCYT= param(45); 
	global_par_Keqcact=1.0;
	global_par_Vcpt2=boxVmax(3);        
	global_par_Kmcpt2C16AcylCarMAT=51.0/SF;
	global_par_Kmcpt2C14AcylCarMAT=51.0/SF;
	global_par_Kmcpt2C12AcylCarMAT=51.0/SF;
	global_par_Kmcpt2C10AcylCarMAT=51.0/SF;
	global_par_Kmcpt2C8AcylCarMAT=51.0/SF;
	global_par_Kmcpt2C6AcylCarMAT=51.0/SF;
	global_par_Kmcpt2C4AcylCarMAT=51.0/SF;
	global_par_Kmcpt2CoAMAT=30.0/SF;
	global_par_Kmcpt2C16AcylCoAMAT=38.0/SF;
	global_par_Kmcpt2C14AcylCoAMAT=38.0/SF;
	global_par_Kmcpt2C12AcylCoAMAT=38.0/SF;
	global_par_Kmcpt2C10AcylCoAMAT=38.0/SF;
	global_par_Kmcpt2C8AcylCoAMAT=38.0/SF;
	global_par_Kmcpt2C6AcylCoAMAT=1000.0/SF;
	global_par_Kmcpt2C4AcylCoAMAT=1000000.0/SF;
	global_par_Kmcpt2CarMAT=350.0/SF;
	global_par_Keqcpt2=2.22;
	global_par_Vvlcad=boxVmax(4);  
	global_par_KmvlcadC16AcylCoAMAT=1.08/SF; 
	global_par_KmvlcadC14AcylCoAMAT=1.08/SF; 
	global_par_KmvlcadC12AcylCoAMAT=1.08/SF; 
	global_par_KmvlcadFAD=0.12/SF;
	global_par_KmvlcadC16EnoylCoAMAT=1.08/SF;
	global_par_KmvlcadC14EnoylCoAMAT=1.08/SF;
	global_par_KmvlcadC12EnoylCoAMAT=1.08/SF;
	global_par_KmvlcadFADH=24.2/SF;
	global_par_Keqvlcad=6.0;
	global_par_Vlcad=boxVmax(5);     
	global_par_KmlcadC16AcylCoAMAT=2.5/SF;
	global_par_KmlcadC14AcylCoAMAT=7.4/SF;
	global_par_KmlcadC12AcylCoAMAT=9.0/SF;
	global_par_KmlcadC10AcylCoAMAT=24.3/SF;
	global_par_KmlcadC8AcylCoAMAT=123.0/SF;
	global_par_KmlcadFAD=0.12/SF;
	global_par_KmlcadC16EnoylCoAMAT=1.08/SF;
	global_par_KmlcadC14EnoylCoAMAT=1.08/SF;
	global_par_KmlcadC12EnoylCoAMAT=1.08/SF;
	global_par_KmlcadC10EnoylCoAMAT=1.08/SF;
	global_par_KmlcadC8EnoylCoAMAT=1.08/SF;
	global_par_KmlcadFADH=24.2/SF;
	global_par_Keqlcad=6.0;
	global_par_Vmcad=boxVmax(6);   
	global_par_KmmcadC12AcylCoAMAT=5.7/SF;
	global_par_KmmcadC10AcylCoAMAT=5.4/SF;
	global_par_KmmcadC8AcylCoAMAT=4.0/SF;
	global_par_KmmcadC6AcylCoAMAT=9.4/SF;
	global_par_KmmcadC4AcylCoAMAT=135.0/SF;
	global_par_KmmcadFAD=0.12/SF;
	global_par_KmmcadC12EnoylCoAMAT=1.08/SF;
	global_par_KmmcadC10EnoylCoAMAT=1.08/SF;
	global_par_KmmcadC8EnoylCoAMAT=1.08/SF;
	global_par_KmmcadC6EnoylCoAMAT=1.08/SF;
	global_par_KmmcadC4EnoylCoAMAT=1.08/SF;
	global_par_KmmcadFADH=24.2/SF;
	global_par_Keqmcad=6.0;
	global_par_Vscad=boxVmax(7);       
	global_par_KmscadC6AcylCoAMAT=285.0/SF;
	global_par_KmscadC4AcylCoAMAT=10.7/SF;
	global_par_KmscadFAD=0.12/SF;
	global_par_KmscadC6EnoylCoAMAT=1.08/SF;
	global_par_KmscadC4EnoylCoAMAT=1.08/SF;
	global_par_KmscadFADH=24.2/SF;
	global_par_Keqscad=6.0;
	global_par_Vcrot=boxVmax(8);     
	global_par_KmcrotC16EnoylCoAMAT=150.0/SF;
	global_par_KmcrotC14EnoylCoAMAT=100.0/SF;
	global_par_KmcrotC12EnoylCoAMAT=25.0/SF;
	global_par_KmcrotC10EnoylCoAMAT=25.0/SF;
	global_par_KmcrotC8EnoylCoAMAT=25.0/SF;
	global_par_KmcrotC6EnoylCoAMAT=25.0/SF;
	global_par_KmcrotC4EnoylCoAMAT=40.0/SF;
	global_par_KmcrotC16HydroxyacylCoAMAT=45.0/SF;
	global_par_KmcrotC14HydroxyacylCoAMAT=45.0/SF;
    global_par_KmcrotC12HydroxyacylCoAMAT=45.0/SF;
	global_par_KmcrotC10HydroxyacylCoAMAT=45.0/SF;
	global_par_KmcrotC8HydroxyacylCoAMAT=45.0/SF;
    global_par_KmcrotC6HydroxyacylCoAMAT=45.0/SF;
	global_par_KmcrotC4HydroxyacylCoAMAT=45.0/SF;
	global_par_KicrotC4AcetoacylCoA=1.6/SF;
	global_par_Keqcrot=3.13;                
	global_par_Vmschad=boxVmax(9);    
    global_par_KmmschadC16HydroxyacylCoAMAT=1.5/SF;
	global_par_KmmschadC14HydroxyacylCoAMAT=1.8/SF;
	global_par_KmmschadC12HydroxyacylCoAMAT=3.7/SF;
	global_par_KmmschadC10HydroxyacylCoAMAT=8.8/SF;
	global_par_KmmschadC8HydroxyacylCoAMAT=16.3/SF;
	global_par_KmmschadC6HydroxyacylCoAMAT=28.6/SF;
	global_par_KmmschadC4HydroxyacylCoAMAT=69.9/SF;
	global_par_KmmschadNADMAT=58.5/SF;
	global_par_KmmschadC16KetoacylCoAMAT=1.4/SF;
	global_par_KmmschadC14KetoacylCoAMAT=1.4/SF;
	global_par_KmmschadC12KetoacylCoAMAT=1.6/SF;
	global_par_KmmschadC10KetoacylCoAMAT=2.3/SF;
	global_par_KmmschadC8KetoacylCoAMAT=4.1/SF;
	global_par_KmmschadC6KetoacylCoAMAT=5.8/SF;
	global_par_KmmschadC4AcetoacylCoAMAT=16.9/SF;
	global_par_KmmschadNADHMAT=5.4/SF;
	global_par_Keqmschad=2.17E-4;       
	global_par_Vmckat=boxVmax(10);
    global_par_KmmckatC16KetoacylCoAMAT=1.1/SF;
	global_par_KmmckatC14KetoacylCoAMAT=1.2/SF;
	global_par_KmmckatC12KetoacylCoAMAT=1.3/SF;
	global_par_KmmckatC10KetoacylCoAMAT=2.1/SF;
	global_par_KmmckatC8KetoacylCoAMAT=3.2/SF;
	global_par_KmmckatC6KetoacylCoAMAT=6.7/SF;
	global_par_KmmckatC4AcetoacylCoAMAT=12.4/SF;
	global_par_KmmckatCoAMAT=26.6/SF;
	global_par_KmmckatC14AcylCoAMAT=13.83/SF;
	global_par_KmmckatC16AcylCoAMAT=13.83/SF;
	global_par_KmmckatC12AcylCoAMAT=13.83/SF;
	global_par_KmmckatC10AcylCoAMAT=13.83/SF;
	global_par_KmmckatC8AcylCoAMAT=13.83/SF;
	global_par_KmmckatC6AcylCoAMAT=13.83/SF;
	global_par_KmmckatC4AcylCoAMAT=13.83/SF;
	global_par_KmmckatAcetylCoAMAT=30.0/SF;
	global_par_Keqmckat=1051.0;
	global_par_Vmtp=boxVmax(11);   
	global_par_KmmtpC16EnoylCoAMAT=25.0/SF;
	global_par_KmmtpC14EnoylCoAMAT=25.0/SF;
	global_par_KmmtpC12EnoylCoAMAT=25.0/SF;
	global_par_KmmtpC10EnoylCoAMAT=25.0/SF;
	global_par_KmmtpC8EnoylCoAMAT=25.0/SF;
	global_par_KmmtpNADMAT=60.0/SF;
	global_par_KmmtpCoAMAT=30.0/SF;
	global_par_KmmtpC14AcylCoAMAT=13.83/SF;
	global_par_KmmtpC16AcylCoAMAT=13.83/SF;
	global_par_KmmtpC12AcylCoAMAT=13.83/SF;
	global_par_KmmtpC10AcylCoAMAT=13.83/SF;
	global_par_KmmtpC8AcylCoAMAT=13.83/SF;
	global_par_KmmtpC6AcylCoAMAT=13.83/SF;
	global_par_KmmtpNADHMAT=50.0/SF;
	global_par_KmmtpAcetylCoAMAT=30.0/SF;
	global_par_Keqmtp=0.71;
    
    
% assignmentRule: variable = CoAMAT
	% CoAMAT=const_species_CoAMATt-(C16CoA_m+C16EnoylCoA_m+C16OHCoA_m+C16KetoCoA_m+C14CoA_m+C14EnoylCoA_m+C14OHCoA_m+C14KetoCoA_m+C12CoA_m+C12EnoylCoA_m+C12OHCoA_m+C12KetoCoA_m+C10CoA_m+C10EnoylCoA_m+C10OHCoA_m+C10KetoCoA_m+C8CoA_m+C8EnoylCoA_m+C8OHCoA_m+C8KetoCoA_m+C6CoA_m+C6EnoylCoA_m+C6OHCoA_m+C6KetoCoA_m+C4CoA_m+C4EnoylCoA_m+C4OHCoA_m+C4KetoCoA_m+AcetylCoAMAT);

    % assignmentRule: variable = C16AcylCoACYT
	%C16AcylCoACYT=(26.8/SF);%*2.71828182845905.^((-0.18/SF3).*t);
    
%% Reactions
	reaction_vcpt1C16_Keqcpt1=0.45;
    reaction_vcpt1C16_Kicpt1MalCoACYT=9.1/SF;
    reaction_vcpt1C16_Kmcpt1C16AcylCarCYT=136.0/SF;
    reaction_vcpt1C16_Kmcpt1C16AcylCoACYT=13.8/SF;
    reaction_vcpt1C16_Kmcpt1CarCYT=125.0/SF;
    reaction_vcpt1C16_Kmcpt1CoACYT=40.7/SF;
    reaction_vcpt1C16_Vcpt1=boxVmax(12);
    reaction_vcpt1C16_ncpt1=2.4799; %Hill coefficient
    reaction_vcpt1C16_sfcpt1C16=1.0;
    reaction_vcpt1C16=function_CPT1(reaction_vcpt1C16_sfcpt1C16, reaction_vcpt1C16_Vcpt1, reaction_vcpt1C16_Kmcpt1C16AcylCoACYT, reaction_vcpt1C16_Kmcpt1CarCYT, reaction_vcpt1C16_Kmcpt1C16AcylCarCYT, reaction_vcpt1C16_Kmcpt1CoACYT, reaction_vcpt1C16_Kicpt1MalCoACYT, reaction_vcpt1C16_Keqcpt1, C16AcylCoACYT, const_species_CarCYT, C16Carn_cy, const_species_CoACYT, const_species_MalCoACYT, reaction_vcpt1C16_ncpt1);
    reaction_vcactC16_KicactC16AcylCarCYT=param(44); %56.0/SF;
    reaction_vcactC16_KmcactC16AcylCarCYT=15.0/SF;
    reaction_vcactC16_KmcactC16AcylCarMAT=15.0/SF;
    reaction_vcactC16=function_CACT(global_par_Vfcact, global_par_Vrcact, reaction_vcactC16_KmcactC16AcylCarCYT, global_par_KmcactCarMAT, reaction_vcactC16_KmcactC16AcylCarMAT, global_par_KmcactCarCYT, reaction_vcactC16_KicactC16AcylCarCYT, global_par_KicactCarCYT, global_par_Keqcact, C16Carn_cy, const_species_CarMAT, C16Carn_m, const_species_CarCYT);
    reaction_vcactC14_KicactC14AcylCarCYT=param(44); %56.0/SF;
    reaction_vcactC14_KmcactC14AcylCarCYT=15.0/SF;
    reaction_vcactC14_KmcactC14AcylCarMAT=15.0/SF;
    reaction_vcactC14=function_CACT(global_par_Vfcact, global_par_Vrcact, reaction_vcactC14_KmcactC14AcylCarCYT, global_par_KmcactCarMAT, reaction_vcactC14_KmcactC14AcylCarMAT, global_par_KmcactCarCYT, reaction_vcactC14_KicactC14AcylCarCYT, global_par_KicactCarCYT, global_par_Keqcact, C14Carn_cy, const_species_CarMAT, C14Carn_m, const_species_CarCYT);
    reaction_vcactC12_KicactC12AcylCarCYT=param(44); %56.0/SF;
    reaction_vcactC12_KmcactC12AcylCarCYT=15.0/SF;
    reaction_vcactC12_KmcactC12AcylCarMAT=15.0/SF;
    reaction_vcactC12=function_CACT(global_par_Vfcact, global_par_Vrcact, reaction_vcactC12_KmcactC12AcylCarCYT, global_par_KmcactCarMAT, reaction_vcactC12_KmcactC12AcylCarMAT, global_par_KmcactCarCYT, reaction_vcactC12_KicactC12AcylCarCYT, global_par_KicactCarCYT, global_par_Keqcact, C12Carn_cy, const_species_CarMAT, C12Carn_m, const_species_CarCYT);
    reaction_vcactC10_KicactC10AcylCarCYT=param(44); %56.0/SF;
    reaction_vcactC10_KmcactC10AcylCarCYT=15.0/SF;
    reaction_vcactC10_KmcactC10AcylCarMAT=15.0/SF;
    reaction_vcactC10=function_CACT(global_par_Vfcact, global_par_Vrcact, reaction_vcactC10_KmcactC10AcylCarCYT, global_par_KmcactCarMAT, reaction_vcactC10_KmcactC10AcylCarMAT, global_par_KmcactCarCYT, reaction_vcactC10_KicactC10AcylCarCYT, global_par_KicactCarCYT, global_par_Keqcact, C10Carn_cy, const_species_CarMAT, C10Carn_m, const_species_CarCYT);
    reaction_vcactC8_KicactC8AcylCarCYT=param(44); %56.0/SF;
    reaction_vcactC8_KmcactC8AcylCarCYT=15.0/SF;
    reaction_vcactC8_KmcactC8AcylCarMAT=15.0/SF;
    reaction_vcactC8=function_CACT(global_par_Vfcact, global_par_Vrcact, reaction_vcactC8_KmcactC8AcylCarCYT, global_par_KmcactCarMAT, reaction_vcactC8_KmcactC8AcylCarMAT, global_par_KmcactCarCYT, reaction_vcactC8_KicactC8AcylCarCYT, global_par_KicactCarCYT, global_par_Keqcact, C8Carn_cy, const_species_CarMAT, C8Carn_m, const_species_CarCYT);
    reaction_vcactC6_KicactC6AcylCarCYT=param(44); %56.0/SF;
    reaction_vcactC6_KmcactC6AcylCarCYT=15.0/SF;
    reaction_vcactC6_KmcactC6AcylCarMAT=15.0/SF;
    reaction_vcactC6=function_CACT(global_par_Vfcact, global_par_Vrcact, reaction_vcactC6_KmcactC6AcylCarCYT, global_par_KmcactCarMAT, reaction_vcactC6_KmcactC6AcylCarMAT, global_par_KmcactCarCYT, reaction_vcactC6_KicactC6AcylCarCYT, global_par_KicactCarCYT, global_par_Keqcact, C6Carn_cy, const_species_CarMAT, C6Carn_m, const_species_CarCYT);
    reaction_vcactC4_KicactC4AcylCarCYT=param(44); %56.0/SF;
    reaction_vcactC4_KmcactC4AcylCarCYT=15.0/SF;
    reaction_vcactC4_KmcactC4AcylCarMAT=15.0/SF;
    reaction_vcactC4=function_CACT(global_par_Vfcact, global_par_Vrcact, reaction_vcactC4_KmcactC4AcylCarCYT, global_par_KmcactCarMAT, reaction_vcactC4_KmcactC4AcylCarMAT, global_par_KmcactCarCYT, reaction_vcactC4_KicactC4AcylCarCYT, global_par_KicactCarCYT, global_par_Keqcact, C4Carn_cy, const_species_CarMAT, C4Carn_m, const_species_CarCYT);
    reaction_vcpt2C16_sfcpt2C16=s_f(1);
    reaction_vcpt2C16=function_CPT2(reaction_vcpt2C16_sfcpt2C16, global_par_Vcpt2, global_par_Kmcpt2C16AcylCarMAT, global_par_Kmcpt2C14AcylCarMAT, global_par_Kmcpt2C12AcylCarMAT, global_par_Kmcpt2C10AcylCarMAT, global_par_Kmcpt2C8AcylCarMAT, global_par_Kmcpt2C6AcylCarMAT, global_par_Kmcpt2C4AcylCarMAT, global_par_Kmcpt2CoAMAT, global_par_Kmcpt2C16AcylCoAMAT, global_par_Kmcpt2C14AcylCoAMAT, global_par_Kmcpt2C12AcylCoAMAT, global_par_Kmcpt2C10AcylCoAMAT, global_par_Kmcpt2C8AcylCoAMAT, global_par_Kmcpt2C6AcylCoAMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2CarMAT, global_par_Keqcpt2, C16Carn_m, C14Carn_m, C12Carn_m, C10Carn_m, C8Carn_m, C6Carn_m, C4Carn_m, CoAMAT, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, const_species_CarMAT);
    reaction_vcpt2C14_sfcpt2C14=s_f(2);
    reaction_vcpt2C14=function_CPT2(reaction_vcpt2C14_sfcpt2C14, global_par_Vcpt2, global_par_Kmcpt2C14AcylCarMAT, global_par_Kmcpt2C16AcylCarMAT, global_par_Kmcpt2C12AcylCarMAT, global_par_Kmcpt2C10AcylCarMAT, global_par_Kmcpt2C8AcylCarMAT, global_par_Kmcpt2C6AcylCarMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2CoAMAT, global_par_Kmcpt2C14AcylCoAMAT, global_par_Kmcpt2C16AcylCoAMAT, global_par_Kmcpt2C12AcylCoAMAT, global_par_Kmcpt2C10AcylCoAMAT, global_par_Kmcpt2C8AcylCoAMAT, global_par_Kmcpt2C6AcylCoAMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2CarMAT, global_par_Keqcpt2, C14Carn_m, C16Carn_m, C12Carn_m, C10Carn_m, C8Carn_m, C6Carn_m, C4Carn_m, CoAMAT, C14CoA_m, C16CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, const_species_CarMAT);
    reaction_vcpt2C12_sfcpt2C12=s_f(3); 
    reaction_vcpt2C12=function_CPT2(reaction_vcpt2C12_sfcpt2C12, global_par_Vcpt2, global_par_Kmcpt2C12AcylCarMAT, global_par_Kmcpt2C16AcylCarMAT, global_par_Kmcpt2C14AcylCarMAT, global_par_Kmcpt2C10AcylCarMAT, global_par_Kmcpt2C8AcylCarMAT, global_par_Kmcpt2C6AcylCarMAT, global_par_Kmcpt2C4AcylCarMAT, global_par_Kmcpt2CoAMAT, global_par_Kmcpt2C12AcylCoAMAT, global_par_Kmcpt2C16AcylCoAMAT, global_par_Kmcpt2C14AcylCoAMAT, global_par_Kmcpt2C10AcylCoAMAT, global_par_Kmcpt2C8AcylCoAMAT, global_par_Kmcpt2C6AcylCoAMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2CarMAT, global_par_Keqcpt2, C12Carn_m, C16Carn_m, C14Carn_m, C10Carn_m, C8Carn_m, C6Carn_m, C4Carn_m, CoAMAT, C12CoA_m, C16CoA_m, C14CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, const_species_CarMAT);
    reaction_vcpt2C10_sfcpt2C10=s_f(4);
    reaction_vcpt2C10=function_CPT2(reaction_vcpt2C10_sfcpt2C10, global_par_Vcpt2, global_par_Kmcpt2C10AcylCarMAT, global_par_Kmcpt2C16AcylCarMAT, global_par_Kmcpt2C14AcylCarMAT, global_par_Kmcpt2C12AcylCarMAT, global_par_Kmcpt2C8AcylCarMAT, global_par_Kmcpt2C6AcylCarMAT, global_par_Kmcpt2C4AcylCarMAT, global_par_Kmcpt2CoAMAT, global_par_Kmcpt2C10AcylCoAMAT, global_par_Kmcpt2C16AcylCoAMAT, global_par_Kmcpt2C14AcylCoAMAT, global_par_Kmcpt2C12AcylCoAMAT, global_par_Kmcpt2C8AcylCoAMAT, global_par_Kmcpt2C6AcylCoAMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2CarMAT, global_par_Keqcpt2, C10Carn_m, C16Carn_m, C14Carn_m, C12Carn_m, C8Carn_m, C6Carn_m, C4Carn_m, CoAMAT, C10CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, const_species_CarMAT);
    reaction_vcpt2C8_sfcpt2C8=delta_sf(1); 
    reaction_vcpt2C8=function_CPT2(reaction_vcpt2C8_sfcpt2C8, global_par_Vcpt2, global_par_Kmcpt2C8AcylCarMAT, global_par_Kmcpt2C16AcylCarMAT, global_par_Kmcpt2C14AcylCarMAT, global_par_Kmcpt2C12AcylCarMAT, global_par_Kmcpt2C10AcylCarMAT, global_par_Kmcpt2C6AcylCarMAT, global_par_Kmcpt2C4AcylCarMAT, global_par_Kmcpt2CoAMAT, global_par_Kmcpt2C8AcylCoAMAT, global_par_Kmcpt2C16AcylCoAMAT, global_par_Kmcpt2C14AcylCoAMAT, global_par_Kmcpt2C12AcylCoAMAT, global_par_Kmcpt2C10AcylCoAMAT, global_par_Kmcpt2C6AcylCoAMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2CarMAT, global_par_Keqcpt2, C8Carn_m, C16Carn_m, C14Carn_m, C12Carn_m, C10Carn_m, C6Carn_m, C4Carn_m, CoAMAT, C8CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C6CoA_m, C4CoA_m, const_species_CarMAT);
    reaction_vcpt2C6_sfcpt2C6=delta_sf(2); 
    reaction_vcpt2C6=function_CPT2(reaction_vcpt2C6_sfcpt2C6, global_par_Vcpt2, global_par_Kmcpt2C6AcylCarMAT, global_par_Kmcpt2C16AcylCarMAT, global_par_Kmcpt2C14AcylCarMAT, global_par_Kmcpt2C12AcylCarMAT, global_par_Kmcpt2C10AcylCarMAT, global_par_Kmcpt2C8AcylCarMAT, global_par_Kmcpt2C4AcylCarMAT, global_par_Kmcpt2CoAMAT, global_par_Kmcpt2C6AcylCoAMAT, global_par_Kmcpt2C16AcylCoAMAT, global_par_Kmcpt2C14AcylCoAMAT, global_par_Kmcpt2C12AcylCoAMAT, global_par_Kmcpt2C10AcylCoAMAT, global_par_Kmcpt2C8AcylCoAMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2CarMAT, global_par_Keqcpt2, C6Carn_m, C16Carn_m, C14Carn_m, C12Carn_m, C10Carn_m, C8Carn_m, C4Carn_m, CoAMAT, C6CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C4CoA_m, const_species_CarMAT);
    reaction_vcpt2C4_sfcpt2C4=s_f(7);                  
    reaction_vcpt2C4=function_CPT2(reaction_vcpt2C4_sfcpt2C4, global_par_Vcpt2, global_par_Kmcpt2C4AcylCarMAT, global_par_Kmcpt2C16AcylCarMAT, global_par_Kmcpt2C14AcylCarMAT, global_par_Kmcpt2C12AcylCarMAT, global_par_Kmcpt2C10AcylCarMAT, global_par_Kmcpt2C8AcylCarMAT, global_par_Kmcpt2C6AcylCarMAT, global_par_Kmcpt2CoAMAT, global_par_Kmcpt2C4AcylCoAMAT, global_par_Kmcpt2C16AcylCoAMAT, global_par_Kmcpt2C14AcylCoAMAT, global_par_Kmcpt2C12AcylCoAMAT, global_par_Kmcpt2C10AcylCoAMAT, global_par_Kmcpt2C8AcylCoAMAT, global_par_Kmcpt2C6AcylCoAMAT, global_par_Kmcpt2CarMAT, global_par_Keqcpt2, C4Carn_m, C16Carn_m, C14Carn_m, C12Carn_m, C10Carn_m, C8Carn_m, C6Carn_m, CoAMAT, C4CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, const_species_CarMAT);
    reaction_vvlcadC16_sfvlcadC16=s_f(8);              
    reaction_vvlcadC16=function_VLCAD(reaction_vvlcadC16_sfvlcadC16, global_par_Vvlcad, global_par_KmvlcadC16AcylCoAMAT, global_par_KmvlcadC14AcylCoAMAT, global_par_KmvlcadC12AcylCoAMAT, global_par_KmvlcadFAD, global_par_KmvlcadC16EnoylCoAMAT, global_par_KmvlcadC14EnoylCoAMAT, global_par_KmvlcadC12EnoylCoAMAT, global_par_KmvlcadFADH, global_par_Keqvlcad, C16CoA_m, C14CoA_m, C12CoA_m, const_species_FADtMAT, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, FADH_m);
    reaction_vvlcadC14_sfvlcadC14=s_f(9);              
    reaction_vvlcadC14=function_VLCAD(reaction_vvlcadC14_sfvlcadC14, global_par_Vvlcad, global_par_KmvlcadC14AcylCoAMAT, global_par_KmvlcadC16AcylCoAMAT, global_par_KmvlcadC12AcylCoAMAT, global_par_KmvlcadFAD, global_par_KmvlcadC14EnoylCoAMAT, global_par_KmvlcadC16EnoylCoAMAT, global_par_KmvlcadC12EnoylCoAMAT, global_par_KmvlcadFADH, global_par_Keqvlcad, C14CoA_m, C16CoA_m, C12CoA_m, const_species_FADtMAT, C14EnoylCoA_m, C16EnoylCoA_m, C12EnoylCoA_m, FADH_m);
    reaction_vvlcadC12_sfvlcadC12=s_f(10);             
    reaction_vvlcadC12=function_VLCAD(reaction_vvlcadC12_sfvlcadC12, global_par_Vvlcad, global_par_KmvlcadC12AcylCoAMAT, global_par_KmvlcadC16AcylCoAMAT, global_par_KmvlcadC14AcylCoAMAT, global_par_KmvlcadFAD, global_par_KmvlcadC12EnoylCoAMAT, global_par_KmvlcadC16EnoylCoAMAT, global_par_KmvlcadC14EnoylCoAMAT, global_par_KmvlcadFADH, global_par_Keqvlcad, C12CoA_m, C16CoA_m, C14CoA_m, const_species_FADtMAT, C12EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, FADH_m);
    reaction_vlcadC16_sflcadC16=delta_sf(3);           
    reaction_vlcadC16=function_LCAD(reaction_vlcadC16_sflcadC16, global_par_Vlcad, global_par_KmlcadC16AcylCoAMAT, global_par_KmlcadC14AcylCoAMAT, global_par_KmlcadC12AcylCoAMAT, global_par_KmlcadC10AcylCoAMAT, global_par_KmlcadC8AcylCoAMAT, global_par_KmlcadFAD, global_par_KmlcadC16EnoylCoAMAT, global_par_KmlcadC14EnoylCoAMAT, global_par_KmlcadC12EnoylCoAMAT, global_par_KmlcadC10EnoylCoAMAT, global_par_KmlcadC8EnoylCoAMAT, global_par_KmlcadFADH, global_par_Keqlcad, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, const_species_FADtMAT, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, FADH_m);
    reaction_vlcadC14_sflcadC14=s_f(12);               
    reaction_vlcadC14=function_LCAD(reaction_vlcadC14_sflcadC14, global_par_Vlcad, global_par_KmlcadC14AcylCoAMAT, global_par_KmlcadC16AcylCoAMAT, global_par_KmlcadC12AcylCoAMAT, global_par_KmlcadC10AcylCoAMAT, global_par_KmlcadC8AcylCoAMAT, global_par_KmlcadFAD, global_par_KmlcadC14EnoylCoAMAT, global_par_KmlcadC16EnoylCoAMAT, global_par_KmlcadC12EnoylCoAMAT, global_par_KmlcadC10EnoylCoAMAT, global_par_KmlcadC8EnoylCoAMAT, global_par_KmlcadFADH, global_par_Keqlcad, C14CoA_m, C16CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, const_species_FADtMAT, C14EnoylCoA_m, C16EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, FADH_m);
    reaction_vlcadC12_sflcadC12=s_f(13);              
    reaction_vlcadC12=function_LCAD(reaction_vlcadC12_sflcadC12, global_par_Vlcad, global_par_KmlcadC12AcylCoAMAT, global_par_KmlcadC16AcylCoAMAT, global_par_KmlcadC14AcylCoAMAT, global_par_KmlcadC10AcylCoAMAT, global_par_KmlcadC8AcylCoAMAT, global_par_KmlcadFAD, global_par_KmlcadC12EnoylCoAMAT, global_par_KmlcadC16EnoylCoAMAT, global_par_KmlcadC14EnoylCoAMAT, global_par_KmlcadC10EnoylCoAMAT, global_par_KmlcadC8EnoylCoAMAT, global_par_KmlcadFADH, global_par_Keqlcad, C12CoA_m, C16CoA_m, C14CoA_m, C10CoA_m, C8CoA_m, const_species_FADtMAT, C14EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, FADH_m);
    reaction_vlcadC10_sflcadC10=s_f(14);              
    reaction_vlcadC10=function_LCAD(reaction_vlcadC10_sflcadC10, global_par_Vlcad, global_par_KmlcadC10AcylCoAMAT, global_par_KmlcadC16AcylCoAMAT, global_par_KmlcadC14AcylCoAMAT, global_par_KmlcadC12AcylCoAMAT, global_par_KmlcadC8AcylCoAMAT, global_par_KmlcadFAD, global_par_KmlcadC10EnoylCoAMAT, global_par_KmlcadC16EnoylCoAMAT, global_par_KmlcadC14EnoylCoAMAT, global_par_KmlcadC12EnoylCoAMAT, global_par_KmlcadC8EnoylCoAMAT, global_par_KmlcadFADH, global_par_Keqlcad, C10CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C8CoA_m, const_species_FADtMAT, C10EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C8EnoylCoA_m, FADH_m);
    reaction_vlcadC8_sflcadC8=s_f(15); 				  
	reaction_vlcadC8=function_LCAD(reaction_vlcadC8_sflcadC8, global_par_Vlcad, global_par_KmlcadC8AcylCoAMAT, global_par_KmlcadC16AcylCoAMAT, global_par_KmlcadC14AcylCoAMAT, global_par_KmlcadC12AcylCoAMAT, global_par_KmlcadC10AcylCoAMAT, global_par_KmlcadFAD, global_par_KmlcadC8EnoylCoAMAT, global_par_KmlcadC16EnoylCoAMAT, global_par_KmlcadC14EnoylCoAMAT, global_par_KmlcadC12EnoylCoAMAT, global_par_KmlcadC10EnoylCoAMAT, global_par_KmlcadFADH, global_par_Keqlcad, C8CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, const_species_FADtMAT, C8EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, FADH_m);
    reaction_vmcadC12_sfmcadC12=s_f(16);              
    reaction_vmcadC12=function_MCAD(reaction_vmcadC12_sfmcadC12, global_par_Vmcad, global_par_KmmcadC12AcylCoAMAT, global_par_KmmcadC10AcylCoAMAT, global_par_KmmcadC8AcylCoAMAT, global_par_KmmcadC6AcylCoAMAT, global_par_KmmcadC4AcylCoAMAT, global_par_KmmcadFAD, global_par_KmmcadC12EnoylCoAMAT, global_par_KmmcadC10EnoylCoAMAT, global_par_KmmcadC8EnoylCoAMAT, global_par_KmmcadC6EnoylCoAMAT, global_par_KmmcadC4EnoylCoAMAT, global_par_KmmcadFADH, global_par_Keqmcad, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, const_species_FADtMAT, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, FADH_m);
    reaction_vmcadC10_sfmcadC10=s_f(17);  
    reaction_vmcadC10=function_MCAD(reaction_vmcadC10_sfmcadC10, global_par_Vmcad, global_par_KmmcadC10AcylCoAMAT, global_par_KmmcadC12AcylCoAMAT, global_par_KmmcadC8AcylCoAMAT, global_par_KmmcadC6AcylCoAMAT, global_par_KmmcadC4AcylCoAMAT, global_par_KmmcadFAD, global_par_KmmcadC10EnoylCoAMAT, global_par_KmmcadC12EnoylCoAMAT, global_par_KmmcadC8EnoylCoAMAT, global_par_KmmcadC6EnoylCoAMAT, global_par_KmmcadC4EnoylCoAMAT, global_par_KmmcadFADH, global_par_Keqmcad, C10CoA_m, C12CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, const_species_FADtMAT, C10EnoylCoA_m, C12EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, FADH_m);
    reaction_vmcadC8_sfmcadC8=delta_sf(4); 
    reaction_vmcadC8=function_MCAD(reaction_vmcadC8_sfmcadC8, global_par_Vmcad, global_par_KmmcadC8AcylCoAMAT, global_par_KmmcadC12AcylCoAMAT, global_par_KmmcadC10AcylCoAMAT, global_par_KmmcadC6AcylCoAMAT, global_par_KmmcadC4AcylCoAMAT, global_par_KmmcadFAD, global_par_KmmcadC8EnoylCoAMAT, global_par_KmmcadC12EnoylCoAMAT, global_par_KmmcadC10EnoylCoAMAT, global_par_KmmcadC6EnoylCoAMAT, global_par_KmmcadC4EnoylCoAMAT, global_par_KmmcadFADH, global_par_Keqmcad, C8CoA_m, C12CoA_m, C10CoA_m, C6CoA_m, C4CoA_m, const_species_FADtMAT, C8EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, FADH_m);
    reaction_vmcadC6_sfmcadC6=delta_sf(5); 
    reaction_vmcadC6=function_MCAD(reaction_vmcadC6_sfmcadC6, global_par_Vmcad, global_par_KmmcadC6AcylCoAMAT, global_par_KmmcadC12AcylCoAMAT, global_par_KmmcadC10AcylCoAMAT, global_par_KmmcadC8AcylCoAMAT, global_par_KmmcadC4AcylCoAMAT, global_par_KmmcadFAD, global_par_KmmcadC6EnoylCoAMAT, global_par_KmmcadC12EnoylCoAMAT, global_par_KmmcadC10EnoylCoAMAT, global_par_KmmcadC8EnoylCoAMAT, global_par_KmmcadC4EnoylCoAMAT, global_par_KmmcadFADH, global_par_Keqmcad, C6CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C4CoA_m, const_species_FADtMAT, C6EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C4EnoylCoA_m, FADH_m);
    reaction_vmcadC4_sfmcadC4=s_f(20); 
    reaction_vmcadC4=function_MCAD(reaction_vmcadC4_sfmcadC4, global_par_Vmcad, global_par_KmmcadC4AcylCoAMAT, global_par_KmmcadC12AcylCoAMAT, global_par_KmmcadC10AcylCoAMAT, global_par_KmmcadC8AcylCoAMAT, global_par_KmmcadC6AcylCoAMAT, global_par_KmmcadFAD, global_par_KmmcadC4EnoylCoAMAT, global_par_KmmcadC12EnoylCoAMAT, global_par_KmmcadC10EnoylCoAMAT, global_par_KmmcadC8EnoylCoAMAT, global_par_KmmcadC6EnoylCoAMAT, global_par_KmmcadFADH, global_par_Keqmcad, C4CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, const_species_FADtMAT, C4EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, FADH_m);
    reaction_vscadC6_sfscadC6=s_f(21); 
    reaction_vscadC6=function_SCAD(reaction_vscadC6_sfscadC6, global_par_Vscad, global_par_KmscadC6AcylCoAMAT, global_par_KmscadC4AcylCoAMAT, global_par_KmscadFAD, global_par_KmscadC6EnoylCoAMAT, global_par_KmscadC4EnoylCoAMAT, global_par_KmscadFADH, global_par_Keqscad, C6CoA_m, C4CoA_m, const_species_FADtMAT, C6EnoylCoA_m, C4EnoylCoA_m, FADH_m);
    reaction_vscadC4_sfscadC4=s_f(22); 
    reaction_vscadC4=function_SCAD(reaction_vscadC4_sfscadC4, global_par_Vscad, global_par_KmscadC4AcylCoAMAT, global_par_KmscadC6AcylCoAMAT, global_par_KmscadFAD, global_par_KmscadC4EnoylCoAMAT, global_par_KmscadC6EnoylCoAMAT, global_par_KmscadFADH, global_par_Keqscad, C4CoA_m, C6CoA_m, const_species_FADtMAT, C4EnoylCoA_m, C6EnoylCoA_m, FADH_m);
    reaction_vcrotC16_sfcrotC16=s_f(23); 
    reaction_vcrotC16=function_CROT(reaction_vcrotC16_sfcrotC16, global_par_Vcrot, global_par_KmcrotC16EnoylCoAMAT, global_par_KmcrotC14EnoylCoAMAT, global_par_KmcrotC12EnoylCoAMAT, global_par_KmcrotC10EnoylCoAMAT, global_par_KmcrotC8EnoylCoAMAT, global_par_KmcrotC6EnoylCoAMAT, global_par_KmcrotC4EnoylCoAMAT, global_par_KmcrotC16HydroxyacylCoAMAT, global_par_KmcrotC14HydroxyacylCoAMAT, global_par_KmcrotC12HydroxyacylCoAMAT, global_par_KmcrotC10HydroxyacylCoAMAT, global_par_KmcrotC8HydroxyacylCoAMAT, global_par_KmcrotC6HydroxyacylCoAMAT, global_par_KmcrotC4HydroxyacylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqcrot, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, C4KetoCoA_m);
    reaction_vcrotC14_sfcrotC14=s_f(24); 
    reaction_vcrotC14=function_CROT(reaction_vcrotC14_sfcrotC14, global_par_Vcrot, global_par_KmcrotC14EnoylCoAMAT, global_par_KmcrotC16EnoylCoAMAT, global_par_KmcrotC12EnoylCoAMAT, global_par_KmcrotC10EnoylCoAMAT, global_par_KmcrotC8EnoylCoAMAT, global_par_KmcrotC6EnoylCoAMAT, global_par_KmcrotC4EnoylCoAMAT, global_par_KmcrotC14HydroxyacylCoAMAT, global_par_KmcrotC16HydroxyacylCoAMAT, global_par_KmcrotC12HydroxyacylCoAMAT, global_par_KmcrotC10HydroxyacylCoAMAT, global_par_KmcrotC8HydroxyacylCoAMAT, global_par_KmcrotC6HydroxyacylCoAMAT, global_par_KmcrotC4HydroxyacylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqcrot, C14EnoylCoA_m, C16EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, C14OHCoA_m, C16OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, C4KetoCoA_m);
    reaction_vcrotC12_sfcrotC12=s_f(25); 
    reaction_vcrotC12=function_CROT(reaction_vcrotC12_sfcrotC12, global_par_Vcrot, global_par_KmcrotC12EnoylCoAMAT, global_par_KmcrotC16EnoylCoAMAT, global_par_KmcrotC14EnoylCoAMAT, global_par_KmcrotC10EnoylCoAMAT, global_par_KmcrotC8EnoylCoAMAT, global_par_KmcrotC6EnoylCoAMAT, global_par_KmcrotC4EnoylCoAMAT, global_par_KmcrotC12HydroxyacylCoAMAT, global_par_KmcrotC16HydroxyacylCoAMAT, global_par_KmcrotC14HydroxyacylCoAMAT, global_par_KmcrotC10HydroxyacylCoAMAT, global_par_KmcrotC8HydroxyacylCoAMAT, global_par_KmcrotC6HydroxyacylCoAMAT, global_par_KmcrotC4HydroxyacylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqcrot, C12EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, C12OHCoA_m, C16OHCoA_m, C14OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, C4KetoCoA_m);
    reaction_vcrotC10_sfcrotC10=s_f(26); 
    reaction_vcrotC10=function_CROT(reaction_vcrotC10_sfcrotC10, global_par_Vcrot, global_par_KmcrotC10EnoylCoAMAT, global_par_KmcrotC16EnoylCoAMAT, global_par_KmcrotC14EnoylCoAMAT, global_par_KmcrotC12EnoylCoAMAT, global_par_KmcrotC8EnoylCoAMAT, global_par_KmcrotC6EnoylCoAMAT, global_par_KmcrotC4EnoylCoAMAT, global_par_KmcrotC10HydroxyacylCoAMAT, global_par_KmcrotC16HydroxyacylCoAMAT, global_par_KmcrotC14HydroxyacylCoAMAT, global_par_KmcrotC12HydroxyacylCoAMAT, global_par_KmcrotC8HydroxyacylCoAMAT, global_par_KmcrotC6HydroxyacylCoAMAT, global_par_KmcrotC4HydroxyacylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqcrot, C10EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, C10OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, C4KetoCoA_m);
    reaction_vcrotC8_sfcrotC8=s_f(27); 
    reaction_vcrotC8=function_CROT(reaction_vcrotC8_sfcrotC8, global_par_Vcrot, global_par_KmcrotC8EnoylCoAMAT, global_par_KmcrotC16EnoylCoAMAT, global_par_KmcrotC14EnoylCoAMAT, global_par_KmcrotC12EnoylCoAMAT, global_par_KmcrotC10EnoylCoAMAT, global_par_KmcrotC6EnoylCoAMAT, global_par_KmcrotC4EnoylCoAMAT, global_par_KmcrotC8HydroxyacylCoAMAT, global_par_KmcrotC16HydroxyacylCoAMAT, global_par_KmcrotC14HydroxyacylCoAMAT, global_par_KmcrotC12HydroxyacylCoAMAT, global_par_KmcrotC10HydroxyacylCoAMAT, global_par_KmcrotC6HydroxyacylCoAMAT, global_par_KmcrotC4HydroxyacylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqcrot, C8EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C6EnoylCoA_m, C4EnoylCoA_m, C8OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C6OHCoA_m, C4OHCoA_m, C4KetoCoA_m);
    reaction_vcrotC6_sfcrotC6=s_f(28); 
    reaction_vcrotC6=function_CROT(reaction_vcrotC6_sfcrotC6, global_par_Vcrot, global_par_KmcrotC6EnoylCoAMAT, global_par_KmcrotC16EnoylCoAMAT, global_par_KmcrotC14EnoylCoAMAT, global_par_KmcrotC12EnoylCoAMAT, global_par_KmcrotC10EnoylCoAMAT, global_par_KmcrotC8EnoylCoAMAT, global_par_KmcrotC4EnoylCoAMAT, global_par_KmcrotC6HydroxyacylCoAMAT, global_par_KmcrotC16HydroxyacylCoAMAT, global_par_KmcrotC14HydroxyacylCoAMAT, global_par_KmcrotC12HydroxyacylCoAMAT, global_par_KmcrotC10HydroxyacylCoAMAT, global_par_KmcrotC8HydroxyacylCoAMAT, global_par_KmcrotC4HydroxyacylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqcrot, C6EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C4EnoylCoA_m, C6OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C4OHCoA_m, C4KetoCoA_m);
    reaction_vcrotC4_sfcrotC4=s_f(29); 
    reaction_vcrotC4=function_CROT(reaction_vcrotC4_sfcrotC4, global_par_Vcrot, global_par_KmcrotC4EnoylCoAMAT, global_par_KmcrotC16EnoylCoAMAT, global_par_KmcrotC14EnoylCoAMAT, global_par_KmcrotC12EnoylCoAMAT, global_par_KmcrotC10EnoylCoAMAT, global_par_KmcrotC8EnoylCoAMAT, global_par_KmcrotC6EnoylCoAMAT, global_par_KmcrotC4HydroxyacylCoAMAT, global_par_KmcrotC16HydroxyacylCoAMAT, global_par_KmcrotC14HydroxyacylCoAMAT, global_par_KmcrotC12HydroxyacylCoAMAT, global_par_KmcrotC10HydroxyacylCoAMAT, global_par_KmcrotC8HydroxyacylCoAMAT, global_par_KmcrotC6HydroxyacylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqcrot, C4EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, C6EnoylCoA_m, C4OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4KetoCoA_m);
    reaction_vmschadC16_sfmschadC16=s_f(30); 
    reaction_vmschadC16=function_MSCHAD(reaction_vmschadC16_sfmschadC16, global_par_Vmschad, global_par_KmmschadC16HydroxyacylCoAMAT, global_par_KmmschadC14HydroxyacylCoAMAT, global_par_KmmschadC12HydroxyacylCoAMAT, global_par_KmmschadC10HydroxyacylCoAMAT, global_par_KmmschadC8HydroxyacylCoAMAT, global_par_KmmschadC6HydroxyacylCoAMAT, global_par_KmmschadC4HydroxyacylCoAMAT, global_par_KmmschadNADMAT, global_par_KmmschadC16KetoacylCoAMAT, global_par_KmmschadC14KetoacylCoAMAT, global_par_KmmschadC12KetoacylCoAMAT, global_par_KmmschadC10KetoacylCoAMAT, global_par_KmmschadC8KetoacylCoAMAT, global_par_KmmschadC6KetoacylCoAMAT, global_par_KmmschadC4AcetoacylCoAMAT, global_par_KmmschadNADHMAT, global_par_Keqmschad, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, const_species_NADtMAT, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, NADHm);
    reaction_vmschadC14_sfmschadC14=s_f(31); 
    reaction_vmschadC14=function_MSCHAD(reaction_vmschadC14_sfmschadC14, global_par_Vmschad, global_par_KmmschadC14HydroxyacylCoAMAT, global_par_KmmschadC16HydroxyacylCoAMAT, global_par_KmmschadC12HydroxyacylCoAMAT, global_par_KmmschadC10HydroxyacylCoAMAT, global_par_KmmschadC8HydroxyacylCoAMAT, global_par_KmmschadC6HydroxyacylCoAMAT, global_par_KmmschadC4HydroxyacylCoAMAT, global_par_KmmschadNADMAT, global_par_KmmschadC14KetoacylCoAMAT, global_par_KmmschadC16KetoacylCoAMAT, global_par_KmmschadC12KetoacylCoAMAT, global_par_KmmschadC10KetoacylCoAMAT, global_par_KmmschadC8KetoacylCoAMAT, global_par_KmmschadC6KetoacylCoAMAT, global_par_KmmschadC4AcetoacylCoAMAT, global_par_KmmschadNADHMAT, global_par_Keqmschad, C14OHCoA_m, C16OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, const_species_NADtMAT, C14KetoCoA_m, C16KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, NADHm);
    reaction_vmschadC12_sfmschadC12=s_f(32); 
    reaction_vmschadC12=function_MSCHAD(reaction_vmschadC12_sfmschadC12, global_par_Vmschad, global_par_KmmschadC12HydroxyacylCoAMAT, global_par_KmmschadC16HydroxyacylCoAMAT, global_par_KmmschadC14HydroxyacylCoAMAT, global_par_KmmschadC10HydroxyacylCoAMAT, global_par_KmmschadC8HydroxyacylCoAMAT, global_par_KmmschadC6HydroxyacylCoAMAT, global_par_KmmschadC4HydroxyacylCoAMAT, global_par_KmmschadNADMAT, global_par_KmmschadC12KetoacylCoAMAT, global_par_KmmschadC16KetoacylCoAMAT, global_par_KmmschadC14KetoacylCoAMAT, global_par_KmmschadC10KetoacylCoAMAT, global_par_KmmschadC8KetoacylCoAMAT, global_par_KmmschadC6KetoacylCoAMAT, global_par_KmmschadC4AcetoacylCoAMAT, global_par_KmmschadNADHMAT, global_par_Keqmschad, C12OHCoA_m, C16OHCoA_m, C14OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, const_species_NADtMAT, C12KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, NADHm);
    reaction_vmschadC10_sfmschadC10=s_f(33); 
    reaction_vmschadC10=function_MSCHAD(reaction_vmschadC10_sfmschadC10, global_par_Vmschad, global_par_KmmschadC10HydroxyacylCoAMAT, global_par_KmmschadC16HydroxyacylCoAMAT, global_par_KmmschadC14HydroxyacylCoAMAT, global_par_KmmschadC12HydroxyacylCoAMAT, global_par_KmmschadC8HydroxyacylCoAMAT, global_par_KmmschadC6HydroxyacylCoAMAT, global_par_KmmschadC4HydroxyacylCoAMAT, global_par_KmmschadNADMAT, global_par_KmmschadC10KetoacylCoAMAT, global_par_KmmschadC16KetoacylCoAMAT, global_par_KmmschadC14KetoacylCoAMAT, global_par_KmmschadC12KetoacylCoAMAT, global_par_KmmschadC8KetoacylCoAMAT, global_par_KmmschadC6KetoacylCoAMAT, global_par_KmmschadC4AcetoacylCoAMAT, global_par_KmmschadNADHMAT, global_par_Keqmschad, C10OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C8OHCoA_m, C6OHCoA_m, C4OHCoA_m, const_species_NADtMAT, C10KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, NADHm);
    reaction_vmschadC8_sfmschadC8=s_f(34); 
    reaction_vmschadC8=function_MSCHAD(reaction_vmschadC8_sfmschadC8, global_par_Vmschad, global_par_KmmschadC8HydroxyacylCoAMAT, global_par_KmmschadC16HydroxyacylCoAMAT, global_par_KmmschadC14HydroxyacylCoAMAT, global_par_KmmschadC12HydroxyacylCoAMAT, global_par_KmmschadC10HydroxyacylCoAMAT, global_par_KmmschadC6HydroxyacylCoAMAT, global_par_KmmschadC4HydroxyacylCoAMAT, global_par_KmmschadNADMAT, global_par_KmmschadC8KetoacylCoAMAT, global_par_KmmschadC16KetoacylCoAMAT, global_par_KmmschadC14KetoacylCoAMAT, global_par_KmmschadC12KetoacylCoAMAT, global_par_KmmschadC10KetoacylCoAMAT, global_par_KmmschadC6KetoacylCoAMAT, global_par_KmmschadC4AcetoacylCoAMAT, global_par_KmmschadNADHMAT, global_par_Keqmschad, C8OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C6OHCoA_m, C4OHCoA_m, const_species_NADtMAT, C8KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, NADHm);
    reaction_vmschadC6_sfmschadC6=s_f(35); 
    reaction_vmschadC6=function_MSCHAD(reaction_vmschadC6_sfmschadC6, global_par_Vmschad, global_par_KmmschadC6HydroxyacylCoAMAT, global_par_KmmschadC16HydroxyacylCoAMAT, global_par_KmmschadC14HydroxyacylCoAMAT, global_par_KmmschadC12HydroxyacylCoAMAT, global_par_KmmschadC10HydroxyacylCoAMAT, global_par_KmmschadC8HydroxyacylCoAMAT, global_par_KmmschadC4HydroxyacylCoAMAT, global_par_KmmschadNADMAT, global_par_KmmschadC6KetoacylCoAMAT, global_par_KmmschadC16KetoacylCoAMAT, global_par_KmmschadC14KetoacylCoAMAT, global_par_KmmschadC12KetoacylCoAMAT, global_par_KmmschadC10KetoacylCoAMAT, global_par_KmmschadC8KetoacylCoAMAT, global_par_KmmschadC4AcetoacylCoAMAT, global_par_KmmschadNADHMAT, global_par_Keqmschad, C6OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C4OHCoA_m, const_species_NADtMAT, C6KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C4KetoCoA_m, NADHm);
    reaction_vmschadC4_sfmschadC4=s_f(36); 
    reaction_vmschadC4=function_MSCHAD(reaction_vmschadC4_sfmschadC4, global_par_Vmschad, global_par_KmmschadC4HydroxyacylCoAMAT, global_par_KmmschadC16HydroxyacylCoAMAT, global_par_KmmschadC14HydroxyacylCoAMAT, global_par_KmmschadC12HydroxyacylCoAMAT, global_par_KmmschadC10HydroxyacylCoAMAT, global_par_KmmschadC8HydroxyacylCoAMAT, global_par_KmmschadC6HydroxyacylCoAMAT, global_par_KmmschadNADMAT, global_par_KmmschadC4AcetoacylCoAMAT, global_par_KmmschadC16KetoacylCoAMAT, global_par_KmmschadC14KetoacylCoAMAT, global_par_KmmschadC12KetoacylCoAMAT, global_par_KmmschadC10KetoacylCoAMAT, global_par_KmmschadC8KetoacylCoAMAT, global_par_KmmschadC6KetoacylCoAMAT, global_par_KmmschadNADHMAT, global_par_Keqmschad, C4OHCoA_m, C16OHCoA_m, C14OHCoA_m, C12OHCoA_m, C10OHCoA_m, C8OHCoA_m, C6OHCoA_m, const_species_NADtMAT, C4KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, NADHm);
    reaction_vmckatC16_sfmckatC16=s_f(37); 
    reaction_vmckatC16=function_MCKATA(reaction_vmckatC16_sfmckatC16, global_par_Vmckat, global_par_KmmckatC16KetoacylCoAMAT, global_par_KmmckatC14KetoacylCoAMAT, global_par_KmmckatC12KetoacylCoAMAT, global_par_KmmckatC10KetoacylCoAMAT, global_par_KmmckatC8KetoacylCoAMAT, global_par_KmmckatC6KetoacylCoAMAT, global_par_KmmckatC4AcetoacylCoAMAT, global_par_KmmckatCoAMAT, global_par_KmmckatC14AcylCoAMAT, global_par_KmmckatC16AcylCoAMAT, global_par_KmmckatC12AcylCoAMAT, global_par_KmmckatC10AcylCoAMAT, global_par_KmmckatC8AcylCoAMAT, global_par_KmmckatC6AcylCoAMAT, global_par_KmmckatC4AcylCoAMAT, global_par_KmmckatAcetylCoAMAT, global_par_Keqmckat, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, CoAMAT, C14CoA_m, C16CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, AcetylCoAMAT);
    reaction_vmckatC14_sfmckatC14=s_f(38); 
    reaction_vmckatC14=function_MCKATA(reaction_vmckatC14_sfmckatC14, global_par_Vmckat, global_par_KmmckatC14KetoacylCoAMAT, global_par_KmmckatC16KetoacylCoAMAT, global_par_KmmckatC12KetoacylCoAMAT, global_par_KmmckatC10KetoacylCoAMAT, global_par_KmmckatC8KetoacylCoAMAT, global_par_KmmckatC6KetoacylCoAMAT, global_par_KmmckatC4AcetoacylCoAMAT, global_par_KmmckatCoAMAT, global_par_KmmckatC12AcylCoAMAT, global_par_KmmckatC16AcylCoAMAT, global_par_KmmckatC14AcylCoAMAT, global_par_KmmckatC10AcylCoAMAT, global_par_KmmckatC8AcylCoAMAT, global_par_KmmckatC6AcylCoAMAT, global_par_KmmckatC4AcylCoAMAT, global_par_KmmckatAcetylCoAMAT, global_par_Keqmckat, C14KetoCoA_m, C16KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, CoAMAT, C12CoA_m, C16CoA_m, C14CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, AcetylCoAMAT);
    reaction_vmckatC12_sfmckatC12=s_f(39);  
    reaction_vmckatC12=function_MCKATA(reaction_vmckatC12_sfmckatC12, global_par_Vmckat, global_par_KmmckatC12KetoacylCoAMAT, global_par_KmmckatC16KetoacylCoAMAT, global_par_KmmckatC14KetoacylCoAMAT, global_par_KmmckatC10KetoacylCoAMAT, global_par_KmmckatC8KetoacylCoAMAT, global_par_KmmckatC6KetoacylCoAMAT, global_par_KmmckatC4AcetoacylCoAMAT, global_par_KmmckatCoAMAT, global_par_KmmckatC10AcylCoAMAT, global_par_KmmckatC16AcylCoAMAT, global_par_KmmckatC14AcylCoAMAT, global_par_KmmckatC12AcylCoAMAT, global_par_KmmckatC8AcylCoAMAT, global_par_KmmckatC6AcylCoAMAT, global_par_KmmckatC4AcylCoAMAT, global_par_KmmckatAcetylCoAMAT, global_par_Keqmckat, C12KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, CoAMAT, C10CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C8CoA_m, C6CoA_m, C4CoA_m, AcetylCoAMAT);
    reaction_vmckatC10_sfmckatC10=s_f(40); 
    reaction_vmckatC10=function_MCKATA(reaction_vmckatC10_sfmckatC10, global_par_Vmckat, global_par_KmmckatC10KetoacylCoAMAT, global_par_KmmckatC16KetoacylCoAMAT, global_par_KmmckatC14KetoacylCoAMAT, global_par_KmmckatC12KetoacylCoAMAT, global_par_KmmckatC8KetoacylCoAMAT, global_par_KmmckatC6KetoacylCoAMAT, global_par_KmmckatC4AcetoacylCoAMAT, global_par_KmmckatCoAMAT, global_par_KmmckatC8AcylCoAMAT, global_par_KmmckatC16AcylCoAMAT, global_par_KmmckatC14AcylCoAMAT, global_par_KmmckatC12AcylCoAMAT, global_par_KmmckatC10AcylCoAMAT, global_par_KmmckatC6AcylCoAMAT, global_par_KmmckatC4AcylCoAMAT, global_par_KmmckatAcetylCoAMAT, global_par_Keqmckat, C10KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, CoAMAT, C8CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C6CoA_m, C4CoA_m, AcetylCoAMAT);
    reaction_vmckatC8_sfmckatC8=s_f(41); 
    reaction_vmckatC8=function_MCKATA(reaction_vmckatC8_sfmckatC8, global_par_Vmckat, global_par_KmmckatC8KetoacylCoAMAT, global_par_KmmckatC16KetoacylCoAMAT, global_par_KmmckatC14KetoacylCoAMAT, global_par_KmmckatC12KetoacylCoAMAT, global_par_KmmckatC10KetoacylCoAMAT, global_par_KmmckatC6KetoacylCoAMAT, global_par_KmmckatC4AcetoacylCoAMAT, global_par_KmmckatCoAMAT, global_par_KmmckatC6AcylCoAMAT, global_par_KmmckatC16AcylCoAMAT, global_par_KmmckatC14AcylCoAMAT, global_par_KmmckatC12AcylCoAMAT, global_par_KmmckatC10AcylCoAMAT, global_par_KmmckatC8AcylCoAMAT, global_par_KmmckatC4AcylCoAMAT, global_par_KmmckatAcetylCoAMAT, global_par_Keqmckat, C8KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C6KetoCoA_m, C4KetoCoA_m, CoAMAT, C6CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C4CoA_m, AcetylCoAMAT);
    reaction_vmckatC6_sfmckatC6=s_f(42); 
    reaction_vmckatC6=function_MCKATA(reaction_vmckatC6_sfmckatC6, global_par_Vmckat, global_par_KmmckatC6KetoacylCoAMAT, global_par_KmmckatC16KetoacylCoAMAT, global_par_KmmckatC14KetoacylCoAMAT, global_par_KmmckatC12KetoacylCoAMAT, global_par_KmmckatC10KetoacylCoAMAT, global_par_KmmckatC8KetoacylCoAMAT, global_par_KmmckatC4AcetoacylCoAMAT, global_par_KmmckatCoAMAT, global_par_KmmckatC4AcylCoAMAT, global_par_KmmckatC16AcylCoAMAT, global_par_KmmckatC14AcylCoAMAT, global_par_KmmckatC12AcylCoAMAT, global_par_KmmckatC10AcylCoAMAT, global_par_KmmckatC8AcylCoAMAT, global_par_KmmckatC6AcylCoAMAT, global_par_KmmckatAcetylCoAMAT, global_par_Keqmckat, C6KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C4KetoCoA_m, CoAMAT, C4CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, AcetylCoAMAT);
    reaction_vmckatC4_sfmckatC4=s_f(43); 
    reaction_vmckatC4=function_MCKATB(reaction_vmckatC4_sfmckatC4, global_par_Vmckat, global_par_KmmckatC4AcetoacylCoAMAT, global_par_KmmckatC16KetoacylCoAMAT, global_par_KmmckatC14KetoacylCoAMAT, global_par_KmmckatC12KetoacylCoAMAT, global_par_KmmckatC10KetoacylCoAMAT, global_par_KmmckatC8KetoacylCoAMAT, global_par_KmmckatC6KetoacylCoAMAT, global_par_KmmckatCoAMAT, global_par_KmmckatC4AcylCoAMAT, global_par_KmmckatC16AcylCoAMAT, global_par_KmmckatC14AcylCoAMAT, global_par_KmmckatC12AcylCoAMAT, global_par_KmmckatC10AcylCoAMAT, global_par_KmmckatC8AcylCoAMAT, global_par_KmmckatC6AcylCoAMAT, global_par_KmmckatAcetylCoAMAT, global_par_Keqmckat, C4KetoCoA_m, C16KetoCoA_m, C14KetoCoA_m, C12KetoCoA_m, C10KetoCoA_m, C8KetoCoA_m, C6KetoCoA_m, CoAMAT, C4CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, AcetylCoAMAT);
    reaction_vmtpC16_sfmtpC16=s_f(44); 
    reaction_vmtpC16=function_MTP(reaction_vmtpC16_sfmtpC16, global_par_Vmtp, global_par_KmmtpC16EnoylCoAMAT, global_par_KmmtpC14EnoylCoAMAT, global_par_KmmtpC12EnoylCoAMAT, global_par_KmmtpC10EnoylCoAMAT, global_par_KmmtpC8EnoylCoAMAT, global_par_KmmtpNADMAT, global_par_KmmtpCoAMAT, global_par_KmmtpC14AcylCoAMAT, global_par_KmmtpC16AcylCoAMAT, global_par_KmmtpC12AcylCoAMAT, global_par_KmmtpC10AcylCoAMAT, global_par_KmmtpC8AcylCoAMAT, global_par_KmmtpC6AcylCoAMAT, global_par_KmmtpNADHMAT, global_par_KmmtpAcetylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqmtp, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, const_species_NADtMAT, CoAMAT, C14CoA_m, C16CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, NADHm, AcetylCoAMAT, C4KetoCoA_m);
    reaction_vmtpC14_sfmtpC14=s_f(45); 
    reaction_vmtpC14=function_MTP(reaction_vmtpC14_sfmtpC14, global_par_Vmtp, global_par_KmmtpC14EnoylCoAMAT, global_par_KmmtpC16EnoylCoAMAT, global_par_KmmtpC12EnoylCoAMAT, global_par_KmmtpC10EnoylCoAMAT, global_par_KmmtpC8EnoylCoAMAT, global_par_KmmtpNADMAT, global_par_KmmtpCoAMAT, global_par_KmmtpC12AcylCoAMAT, global_par_KmmtpC16AcylCoAMAT, global_par_KmmtpC14AcylCoAMAT, global_par_KmmtpC10AcylCoAMAT, global_par_KmmtpC8AcylCoAMAT, global_par_KmmtpC6AcylCoAMAT, global_par_KmmtpNADHMAT, global_par_KmmtpAcetylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqmtp, C14EnoylCoA_m, C16EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, const_species_NADtMAT, CoAMAT, C12CoA_m, C16CoA_m, C14CoA_m, C10CoA_m, C8CoA_m, C6CoA_m, NADHm, AcetylCoAMAT, C4KetoCoA_m);
    reaction_vmtpC12_sfmtpC12=s_f(46); 
    reaction_vmtpC12=function_MTP(reaction_vmtpC12_sfmtpC12, global_par_Vmtp, global_par_KmmtpC12EnoylCoAMAT, global_par_KmmtpC16EnoylCoAMAT, global_par_KmmtpC14EnoylCoAMAT, global_par_KmmtpC10EnoylCoAMAT, global_par_KmmtpC8EnoylCoAMAT, global_par_KmmtpNADMAT, global_par_KmmtpCoAMAT, global_par_KmmtpC10AcylCoAMAT, global_par_KmmtpC16AcylCoAMAT, global_par_KmmtpC14AcylCoAMAT, global_par_KmmtpC12AcylCoAMAT, global_par_KmmtpC8AcylCoAMAT, global_par_KmmtpC6AcylCoAMAT, global_par_KmmtpNADHMAT, global_par_KmmtpAcetylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqmtp, C12EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C10EnoylCoA_m, C8EnoylCoA_m, const_species_NADtMAT, CoAMAT, C10CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C8CoA_m, C6CoA_m, NADHm, AcetylCoAMAT, C4KetoCoA_m);
    reaction_vmtpC10_sfmtpC10=s_f(47); 
    reaction_vmtpC10=function_MTP(reaction_vmtpC10_sfmtpC10, global_par_Vmtp, global_par_KmmtpC10EnoylCoAMAT, global_par_KmmtpC16EnoylCoAMAT, global_par_KmmtpC14EnoylCoAMAT, global_par_KmmtpC12EnoylCoAMAT, global_par_KmmtpC8EnoylCoAMAT, global_par_KmmtpNADMAT, global_par_KmmtpCoAMAT, global_par_KmmtpC8AcylCoAMAT, global_par_KmmtpC16AcylCoAMAT, global_par_KmmtpC14AcylCoAMAT, global_par_KmmtpC12AcylCoAMAT, global_par_KmmtpC10AcylCoAMAT, global_par_KmmtpC6AcylCoAMAT, global_par_KmmtpNADHMAT, global_par_KmmtpAcetylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqmtp, C10EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C8EnoylCoA_m, const_species_NADtMAT, CoAMAT, C8CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C6CoA_m, NADHm, AcetylCoAMAT, C4KetoCoA_m);
    reaction_vmtpC8_sfmtpC8=s_f(48); 
    reaction_vmtpC8=function_MTP(reaction_vmtpC8_sfmtpC8, global_par_Vmtp, global_par_KmmtpC8EnoylCoAMAT, global_par_KmmtpC16EnoylCoAMAT, global_par_KmmtpC14EnoylCoAMAT, global_par_KmmtpC12EnoylCoAMAT, global_par_KmmtpC10EnoylCoAMAT, global_par_KmmtpNADMAT, global_par_KmmtpCoAMAT, global_par_KmmtpC6AcylCoAMAT, global_par_KmmtpC16AcylCoAMAT, global_par_KmmtpC14AcylCoAMAT, global_par_KmmtpC12AcylCoAMAT, global_par_KmmtpC10AcylCoAMAT, global_par_KmmtpC8AcylCoAMAT, global_par_KmmtpNADHMAT, global_par_KmmtpAcetylCoAMAT, global_par_KicrotC4AcetoacylCoA, global_par_Keqmtp, C8EnoylCoA_m, C16EnoylCoA_m, C14EnoylCoA_m, C12EnoylCoA_m, C10EnoylCoA_m, const_species_NADtMAT, CoAMAT, C6CoA_m, C16CoA_m, C14CoA_m, C12CoA_m, C10CoA_m, C8CoA_m, NADHm, AcetylCoAMAT, C4KetoCoA_m);
    reaction_vacesink_K1acesink=30.0/SF;
    reaction_vacesink_Ksacesink=6000000.0/SF/SF3;
    reaction_vacesink=function_RES(reaction_vacesink_Ksacesink, AcetylCoAMAT, reaction_vacesink_K1acesink);
    reaction_vfadhsink_K1fadhsink=0.46/SF;
    reaction_vfadhsink_Ksfadhsink=6000000.0/SF/SF3;
    reaction_vfadhsink=function_RES(reaction_vfadhsink_Ksfadhsink, FADH_m, reaction_vfadhsink_K1fadhsink);
    reaction_vnadhsink_K1nadhsink= 1.6584E-6; 
    reaction_vnadhsink_Ksnadhsink=6000000.0/SF/SF3;
    reaction_vnadhsink=function_RES(reaction_vnadhsink_Ksnadhsink, NADHm, reaction_vnadhsink_K1nadhsink);
    KsCoASH = 6000000.0/SF/SF3;
    K1CoASH = 9.4686e-6;
    reaction_CoASHsink = function_RES(KsCoASH,CoAMAT,K1CoASH);


%% Beta oxidation species

	f(iC16Carn_cy) = (1/(compartment_VCYT))*(( 1.0 * reaction_vcpt1C16) + (-1.0 * reaction_vcactC16));
	f(iC16Carn_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcactC16) + (-1.0 * reaction_vcpt2C16));
	f(iC16CoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcpt2C16) + (-1.0 * reaction_vvlcadC16) + (-1.0 * reaction_vlcadC16));
	f(iC16EnoylCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vvlcadC16) + ( 1.0 * reaction_vlcadC16) + (-1.0 * reaction_vcrotC16) + (-1.0 * reaction_vmtpC16));
	f(iC16OHCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcrotC16) + (-1.0 * reaction_vmschadC16));
	f(iC16KetoCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC16) + (-1.0 * reaction_vmckatC16));
	f(iC14Carn_cy)= (1/(compartment_VCYT))*((-1.0 * reaction_vcactC14));
	f(iC14Carn_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcactC14) + (-1.0 * reaction_vcpt2C14));
    f(iC14CoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcpt2C14) + (-1.0 * reaction_vvlcadC14) + (-1.0 * reaction_vlcadC14) + ( 1.0 * reaction_vmckatC16) + ( 1.0 * reaction_vmtpC16));
	f(iC14EnoylCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vvlcadC14) + ( 1.0 * reaction_vlcadC14) + (-1.0 * reaction_vcrotC14) + (-1.0 * reaction_vmtpC14));
	f(iC14OHCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcrotC14) + (-1.0 * reaction_vmschadC14));
	f(iC14KetoCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC14) + (-1.0 * reaction_vmckatC14));
	f(iC12Carn_cy) = (1/(compartment_VCYT))*((-1.0 * reaction_vcactC12));
	f(iC12Carn_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcactC12) + (-1.0 * reaction_vcpt2C12));
	f(iC12CoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcpt2C12) + (-1.0 * reaction_vvlcadC12) + (-1.0 * reaction_vlcadC12) + (-1.0 * reaction_vmcadC12) + ( 1.0 * reaction_vmckatC14) + ( 1.0 * reaction_vmtpC14));
	f(iC12EnoylCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vvlcadC12) + ( 1.0 * reaction_vlcadC12) + ( 1.0 * reaction_vmcadC12) + (-1.0 * reaction_vcrotC12) + (-1.0 * reaction_vmtpC12));
	f(iC12OHCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcrotC12) + (-1.0 * reaction_vmschadC12)); %
	f(iC12KetoCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC12) + (-1.0 * reaction_vmckatC12));%
	f(iC10Carn_cy) = (1/(compartment_VCYT))*((-1.0 * reaction_vcactC10));
	f(iC10Carn_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcactC10) + (-1.0 * reaction_vcpt2C10));
	f(iC10CoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcpt2C10) + (-1.0 * reaction_vlcadC10) + (-1.0 * reaction_vmcadC10) + ( 1.0 * reaction_vmckatC12) + ( 1.0 * reaction_vmtpC12));
	f(iC10EnoylCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vlcadC10) + ( 1.0 * reaction_vmcadC10) + (-1.0 * reaction_vcrotC10) + (-1.0 * reaction_vmtpC10));%
	f(iC10OHCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcrotC10) + (-1.0 * reaction_vmschadC10));%
	f(iC10KetoCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC10) + (-1.0 * reaction_vmckatC10));%
	f(iC8Carn_cy) = (1/(compartment_VCYT))*((-1.0 * reaction_vcactC8));
	f(iC8Carn_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcactC8) + (-1.0 * reaction_vcpt2C8));
	f(iC8CoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcpt2C8) + (-1.0 * reaction_vlcadC8) + (-1.0 * reaction_vmcadC8) + ( 1.0 * reaction_vmckatC10) + ( 1.0 * reaction_vmtpC10));
	f(iC8EnoylCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vlcadC8) + ( 1.0 * reaction_vmcadC8) + (-1.0 * reaction_vcrotC8) + (-1.0 * reaction_vmtpC8));%
	f(iC8OHCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcrotC8) + (-1.0 * reaction_vmschadC8));%
	f(iC8KetoCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC8) + (-1.0 * reaction_vmckatC8));%
	f(iC6Carn_cy) = (1/(compartment_VCYT))*((-1.0 * reaction_vcactC6));
	f(iC6Carn_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcactC6) + (-1.0 * reaction_vcpt2C6));
	f(iC6CoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcpt2C6) + (-1.0 * reaction_vmcadC6) + (-1.0 * reaction_vscadC6) + ( 1.0 * reaction_vmckatC8) + ( 1.0 * reaction_vmtpC8));
	f(iC6EnoylCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmcadC6) + ( 1.0 * reaction_vscadC6) + (-1.0 * reaction_vcrotC6));%
	f(iC6OHCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcrotC6) + (-1.0 * reaction_vmschadC6));%
	f(iC6KetoCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC6) + (-1.0 * reaction_vmckatC6));%
	f(iC4Carn_cy) = (1/(compartment_VCYT))*((-1.0 * reaction_vcactC4));
	f(iC4Carn_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcactC4) + (-1.0 * reaction_vcpt2C4));
	f(iC4CoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcpt2C4) + (-1.0 * reaction_vmcadC4) + (-1.0 * reaction_vscadC4) + ( 1.0 * reaction_vmckatC6));%
	f(iC4EnoylCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmcadC4) + ( 1.0 * reaction_vscadC4) + (-1.0 * reaction_vcrotC4));%
	f(iC4OHCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vcrotC4) + (-1.0 * reaction_vmschadC4));
	f(iC4KetoCoA_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC4) + (-1.0 * reaction_vmckatC4));
	%f(iAcetylCoAMAT) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmckatC16) + ( 1.0 * reaction_vmckatC14) + ( 1.0 * reaction_vmckatC12) + ( 1.0 * reaction_vmckatC10) + ( 1.0 * reaction_vmckatC8) + ( 1.0 * reaction_vmckatC6) + ( 2.0 * reaction_vmckatC4) + ( 1.0 * reaction_vmtpC16) + ( 1.0 * reaction_vmtpC14) + ( 1.0 * reaction_vmtpC12) + ( 1.0 * reaction_vmtpC10) + ( 1.0 * reaction_vmtpC8) + (-1.0 * reaction_vacesink));
	f(iFADH_m) = (1/(compartment_VMAT))*(( 1.0 * reaction_vvlcadC16) + ( 1.0 * reaction_vvlcadC14) + ( 1.0 * reaction_vvlcadC12) + ( 1.0 * reaction_vlcadC16) + ( 1.0 * reaction_vlcadC14) + ( 1.0 * reaction_vlcadC12) + ( 1.0 * reaction_vlcadC10) + ( 1.0 * reaction_vlcadC8) + ( 1.0 * reaction_vmcadC12) + ( 1.0 * reaction_vmcadC10) + ( 1.0 * reaction_vmcadC8) + ( 1.0 * reaction_vmcadC6) + ( 1.0 * reaction_vmcadC4) + ( 1.0 * reaction_vscadC6) + ( 1.0 * reaction_vscadC4) + (-1.0 * reaction_vfadhsink));
	%f(iNADHm) = (1/(compartment_VMAT))*(( 1.0 * reaction_vmschadC16) + ( 1.0 * reaction_vmschadC14) + ( 1.0 * reaction_vmschadC12) + ( 1.0 * reaction_vmschadC10) + ( 1.0 * reaction_vmschadC8) + ( 1.0 * reaction_vmschadC6) + ( 1.0 * reaction_vmschadC4) + ( 1.0 * reaction_vmtpC16) + ( 1.0 * reaction_vmtpC14) + ( 1.0 * reaction_vmtpC12) + ( 1.0 * reaction_vmtpC10) + ( 1.0 * reaction_vmtpC8) + (-1.0 * reaction_vnadhsink));
	f(iCoAMAT) = CoAMAT;
    f(iC16AcylCoACYT) = 0;
%% Computing time derivatives of state variables
J_boxaccoa= ((1/compartment_VMAT)* ((reaction_vmckatC16)+(reaction_vmckatC14) + (reaction_vmckatC12) + ( reaction_vmckatC10) + (reaction_vmckatC8) + (reaction_vmckatC6) + ( 2.0 * reaction_vmckatC4) + (  reaction_vmtpC16) + ( reaction_vmtpC14) + (reaction_vmtpC12) + (reaction_vmtpC10) + (reaction_vmtpC8)))/W_x ;   %AcCoA production from beta oxidation
J_boxnadh= ((1/compartment_VMAT)* (( 1.0 * reaction_vmschadC16)+ ( 1.0 * reaction_vmschadC14) + ( 1.0 * reaction_vmschadC12) + ( 1.0 * reaction_vmschadC10) + ( 1.0 * reaction_vmschadC8) + ( 1.0 * reaction_vmschadC6) + ( 1.0 * reaction_vmschadC4) + ( 1.0 * reaction_vmtpC16) + ( 1.0 * reaction_vmtpC14) + ( 1.0 * reaction_vmtpC12) + ( 1.0 * reaction_vmtpC10) + ( 1.0 * reaction_vmtpC8)))/W_x;  %+(-1.0 * reaction_vnadhsink)


%% Computing time derivatives of state variables
% Time derivatives:

% %  (0) cell oxygen:
% if iflag == 1
%   f(iPO2)  = (-J_C4/2)*Rm_cell/(a_3 + CMb*P50/(PO2 + P50)^2); % mol s^{-1} (l cell)^{-1}
% else
%   f(iPO2)  = 0;
% end

% Constant oxygenc
f(iPO2)  = 0;

%  (i) Matrix species and dPsi
f(idPsi)   = (1.48e4)*( 4*J_C1 + 2*J_C3 + 4*J_C4 - n_A*J_F1 - J_ANT - J_Hle + J_ASP_GLU);
f(iATP_x)  = (+J_ndk + J_F1 - J_ANT)/W_x; 
f(iADP_x)  = (-J_ndk - J_F1 + J_ANT)/W_x;
f(iAMP_x)  = (0)/W_x; 
f(iGTP_x)  = (+J_scoas - J_ndk)/W_x; 
f(iGDP_x)  = (-J_scoas + J_ndk)/W_x; 
f(iPI_x)   = (-J_scoas - J_F1 + J_PI1  - J_MAL_PI)/W_x;
f(iNADH_x) = ((+J_pdh + J_isod + J_akgd + J_mdh - J_C1)/W_x)+ J_boxnadh ;
f(iQH2_x)  = (+J_sdh + J_C1 - J_C3)/W_x; 
f(iPYR_x)  = (-J_pdh + J_PYR_H)/W_x; 

f(iACCOA_x) =((-J_cits + J_pdh)/W_x)+ J_boxaccoa; 

f(iCIT_x)  = (+J_cits - J_acon + J_CIT_MAL)/W_x; 
f(iICIT_x) = (+J_acon - J_isod)/W_x; 
f(iAKG_x)  = (+J_isod - J_akgd - J_got + J_AKG_MAL)/W_x; 
f(iSCOA_x) = (+J_akgd - J_scoas)/W_x; 

f(iSUC_x)  = (+J_scoas - J_sdh + J_SUC_MAL)/W_x;
f(iFUM_x)  = (+J_sdh - J_fum )/W_x;  
f(iMAL_x)  = (+J_fum - J_mdh + J_MAL_PI - J_AKG_MAL - J_CIT_MAL - J_SUC_MAL)/W_x; 
f(iOAA_x)  = (-J_cits + J_mdh + J_got)/W_x;
f(iGLU_x)  = (+J_got + J_GLU_H - J_ASP_GLU)/W_x; 
f(iASP_x)  = (-J_got + J_ASP_GLU)/W_x; 
f(iCO2tot_x) = 0;

%  (ii) IM space species
f(iCred_i) = (+2*J_C3 - 2*J_C4)/W_i; 
f(iATP_i) = (J_ATP + J_ANT )/W_i; 
f(iADP_i) = (J_ADP - J_ANT )/W_i; 
f(iAMP_i) = (J_AMP )/W_i;   
f(iPI_i)  = (-J_PI1 + J_PI2 + J_MAL_PI)/W_i; 
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
f(iCARN_i)=J_CARN/W_i;


%  (iii) Buffer species
f(iPYR_c) = (-Rm_buffer*J_PYRt)/W_c; 
f(iCIT_c) = (-Rm_buffer*J_CITt)/W_c;
f(iICIT_c) = (-Rm_buffer*J_ICITt)/W_c;
f(iAKG_c) = (-Rm_buffer*J_AKGt)/W_c;
f(iSUC_c) = (-Rm_buffer*J_SUCt)/W_c;
f(iFUM_c) = (-Rm_buffer*J_FUMt)/W_c;
f(iMAL_c) = (-Rm_buffer*J_MALt)/W_c;
f(iGLU_c) = (-Rm_buffer*J_GLUt)/W_c;
f(iASP_c) = (-Rm_buffer*J_ASPt)/W_c;
f(iH_c)   = 0;
f(iMg_c)  = 0;
f(iK_c)   = 0;
f(iPI_c)  = (-Rm_buffer*J_PI2)/W_c; 
f(iATP_c) = (-Rm_buffer*J_ATP)/W_c;  
f(iADP_c) = (-Rm_buffer*J_ADP)/W_c; 
f(iPCr_c) = 0;
f(iAMP_c) = (-Rm_buffer*J_AMP)/W_c;
f(iCr_c) = 0;
f(iCARN_c)= -(Rm_buffer*J_CARN)/W_c;


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
xx(iPI)     =   x(iPI_x);
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
xx(iFADH2)  =   FADtot/2; 
xx(iFAD)    =   FADtot - xx(iFADH2); 
xx(iCO2tot) =   x(iCO2tot_x); 

dxxdt(1:N_reactant) = 0;
dxxdt(iH)   = f(iH_x);
dxxdt(iATP) = f(iATP_x);
dxxdt(iADP) = f(iADP_x);
dxxdt(iAMP) = f(iAMP_x);
dxxdt(iGTP) = f(iGTP_x);
dxxdt(iGDP) = f(iGDP_x);
dxxdt(iPI)  = f(iPI_x);
dxxdt(iNADH)= f(iNADH_x);
dxxdt(iNAD) = - f(iNADH_x);
dxxdt(iQH2) = f(iQH2_x);
dxxdt(iCOQ) = - f(iQH2_x);
dxxdt(iPYR) = -f(iPYR_x);
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
dxxdt(iK)   = f(iK_x);
dxxdt(iMg)  = f(iMg_x);
dxxdt(iCred)= f(iCred_i);
dxxdt(iCox) = - f(iCred_i);
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
           + 2*J_PI1 + J_Hle - 1*J_KH)/W_x;
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

f(iH_x) = 0; %dH_xdt;
f(iMg_x) = 0; %dMg_xdt;
f(iK_x) = 0; %dK_xdt;
%% Output the time derivative matrix
f = f';
end