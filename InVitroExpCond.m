function InVitroExpCond(fig_flag,StateType,ExpType)

% Copyright: 
% Fan Wu and Daniel. A Beard, May 2007
% Computational Engineering Group
% Biotechnology and Bioengineering Center
% Medical College of Wisconsin

%% Defining global variables% close all;
clear global;

global iH_x idPsi iATP_x iADP_x iAMP_x iGTP_x iGDP_x iPi_x iNADH_x ...
       iQH2_x iOAA_x iACCOA_x iCIT_x iICIT_x iAKG_x iSCOA_x iCOASH_x ...
       iSUC_x iFUM_x iMAL_x iGLU_x iASP_x iK_x iMg_x iO2_x iCO2tot_x ...
       iCred_i iATP_i iADP_i iAMP_i iPi_i iH_i iMg_i iK_i ...
       iATP_c iADP_c iPi_c iH_c iMg_c iK_c ...
       iPYR_x iPYR_i iPYR_c iCIT_i iCIT_c iAKG_i iAKG_c iSUC_i iSUC_c ...
       iMAL_i iMAL_c iASP_i iASP_c iGLU_i iGLU_c iFUM_i iFUM_c ...
       iICIT_i iICIT_c iGLC_c iG6P_c iPCr_c iAMP_c;

define_global1();

%% Setting physical parameters and fixed pools
W_m = 0.664*1.09;   % mitochondrial water space (ml water per ml mito) [VB2002]
W_x = 0.9*W_m;      % Matrix water space, (l mito water) (l mito)^{-1}  (90% of total mito water space is matrix)
W_i = 0.1*W_m;      % IM water space, (l IM water) (l mito)^{-1}  (10% of total mito water space is intermembrane)
rho_m = 3.6697e-6;  % (l mito) (mg protein)^{-1} 
mg_protein=0.025*2; % assume 0.025 mg/mL * 2mL rxn volume %Pierce BCA assay
L_mito=rho_m*mg_protein; %(assumes 2 mL rxn volume)
W_c = (2/1000)/L_mito;  % buffer, (l buffer water) (l mito)^{-1}  assumes 2 mL rxn volume

Ctot   = 2.70e-3;               % M; total cytoC
Qtot   = 1.35e-3;               % M; total Q + QH2
NADtot = 2.97e-3;               % M; total NAD+NADH
% Define Conversion Factor CF to go from nmol/cytA to nmol/mg or nmol/mL
CF= (1/W_c)*1e9*1.8; %  pmol/sec/mL
CF_mg=CF/rho_m/1000; % pmol/sec/mg

%% Loading initial conditions
xo = default_IC_vitro();
% Defining the minimum concentration allowed for components
MinCon = 1e-32; 

%% Loading paramters and experimental data
load param14.mat;

%% Running simulations
indexC = [1 3:63];
options = odeset('MaxStep',5e-2,'NonNegative', indexC);
%% Buffer D: Potassium-MES (105 mM; pH=7.2; KCl 30 mM; KH2PO4 10 mM; MgCL2 5 mM;..
% EGTA 1 mM; BSA 2.5 g/L; ATP 5 mM; Cr 5 mM; PCr 1 mM; CK 20 U/mL; mito 0.025 mg/mL...
% Glutamate/Malate 10/2.5 mM; Pyruvate/Malate 5/2.5 mM; Succinate/rotenone 10/0.005 mM)
% buffer conditions
xo(iH_c)   = 10^(-7.2); 
xo(iMg_c)  = 5e-3; % M 
xo(iK_c)   = 0.150; % M
xo(iH_i)   = xo(iH_c);
xo(iMg_i)  = xo(iMg_c);
xo(iK_i)   = xo(iK_c);
xo(iH_x)   =10^(-7.8);
xo(iATP_c) = 0.005 ; % M
xo(iPCr_c) = 1e-3; % M
xo(iGLC_c)  =10e-3; % 10 mM glucose
xo(iADP_c) = 0.2e-3; % 0.2 mM ADP
xo(iPi_c)  = 20e-3;%  M




%% Figure 2
if fig_flag ==2
%Initialize SS conditions save new initial conditions State 2
valy=ones(1,11) ;
valy=valy*30 ;
valx=[1;1;1;3;3;3;6;6;6;9;12;15;18;21;24;27;30;];
valx=[valx',valy]./1000; %conversion from mM to M

    for i=1:1:2
    xo(iGLU_c) = 10e-6; %5e-3; 10 mM glutamate
    xo(iMAL_c) = 25e-6; %5e-3; 25 mM
    
     xo(iPi_c)  = 20e-3;
     xo(iADP_c) = 0.2e-3;   
J_AtC = .015e-3; %M/sec/mg
[t,x] = ode15s(@IsoMito_dXdT,[0 60],xo,options,param,ExpType,StateType,J_AtC);
xo = x(end,:);

    H_x_PI(i) = x(end,iH_x);
    NADH_x_PI(i) = x(end,iNADH_x);
    Cred_i_PI(i) = x(end,iCred_i);
    dPsi_PI(i)  =x(end,idPsi);       
    ATP_x_PI(i) =x(end,iATP_x);
    ADP_x_PI(i) =x(end,iADP_x);
    ATP_c_PI(i) =x(end,iATP_c);
    ADP_c_PI(i) =x(end,iADP_c);
    PCr_c_PI(i) =x(end,iPCr_c);
    Pi_x_PI(i) =x(end,iPi_x);
    Pi_c_PI(i) =x(end,iPi_c);
  
    t_new(i) =i*t(end);
    xo_con(i,:)=x(end,:);
    t2(i,:)=t(i,:);
   
    t_beg(i)=i*60-(60);
    
    
    
    % Computing oxygen consumption rates
    J = fluxes(x(end,:),param,ExpType,StateType);
    JO2_r(i) = (J(3)/2)*CF; % CF is conversion to pmol/sec/mL; % 1.8957e5; % Unit conversion factor 
    end
 
%% Plotting ready-to-publish curves


      xo(iH_c)   = 10^(-7.2); %org. 7.1xo(iMg_c)  = 5e-3; %M
      xo(iK_c)   = 0.150; %150e-3; M
      xo(iMg_c)  = 5e-3; %M
      xo(iH_i)   = xo(iH_c);
      xo(iMg_i)  = xo(iMg_c);
      xo(iK_i)   = xo(iK_c);
      xo(iH_x)   =10^(-7.8);
      xo(iGLC_c)  =10e-3; % 10 mM glucose
      xo(iADP_c) = 0.2e-3; % 0.2 mM ADP
      xo(iATP_c)= .005 ; %M
      J_AtC = .015e-3; %mM/sec/mg ref. Muoio . (JATP G/M heart 5000 pmol/sec/mg)

    for i=3:1:5
   
    % Adding glutamate and malate
    xo(iGLU_c) = 0.01;   % 10 mM glutamate
    xo(iMAL_c) = 0.0025; % 25 mM
    xo(iPi_c)  = 20e-3;
    xo(iADP_c) = 0.2e-3;    
    [t,x] = ode15s(@IsoMito_dXdT,[0 60],xo,options,param,ExpType,StateType,J_AtC);
    xo1 = x(end,:);
    
    H_x_PI(i) = x(end,iH_x);
    NADH_x_PI(i) = x(end,iNADH_x);
    Cred_i_PI(i) = x(end,iCred_i);
    dPsi_PI(i)  =x(end,idPsi);       
    ATP_x_PI(i) =x(end,iATP_x);
    ADP_x_PI(i) =x(end,iADP_x);
    ATP_c_PI(i) =x(end,iATP_c);
    ADP_c_PI(i) =x(end,iADP_c);
    PCr_c_PI(i) =x(end,iPCr_c);
    Pi_x_PI(i) =x(end,iPi_x);
    Pi_c_PI(i) =x(end,iPi_c);
  
    t_new(i) =i*t(end);
    t_beg(i)=i*60-(60);
    xo_con(i,:)=x(end,:);
    t2(i,:)=t(i,:);
    
    % Computing oxygen consumption rates
    J = fluxes(x(end,:),param,ExpType,StateType);
    JO2_r(i) = (J(3)/2)*CF; %1.8957e5; % Unit conversion factor 
    end   
      
    for i = 6:1:16
      
    % Adding ADP into the buffer
    xo1(iADP_c) = 0.2e-3;    
    PCr_a(i) = valx(i-5);
    xo1(iPCr_c) = PCr_a(i); 
    
    [t,x] = ode15s(@IsoMito_dXdT,[0 60],xo1,options,param,ExpType,StateType,J_AtC);
    xo1 = x(end,:);
    
    H_x_PI(i) = x(end,iH_x);
    NADH_x_PI(i) = x(end,iNADH_x);
    Cred_i_PI(i) = x(end,iCred_i);
    dPsi_PI(i)  =x(end,idPsi);       
    ATP_x_PI(i) =x(end,iATP_x);
    ADP_x_PI(i) =x(end,iADP_x);
    ATP_c_PI(i) =x(end,iATP_c);
    ADP_c_PI(i) =x(end,iADP_c);
    PCr_c_PI(i) =x(end,iPCr_c);
    Pi_x_PI(i) =x(end,iPi_x);
    Pi_c_PI(i) =x(end,iPi_c);
  
    t_new(i) =i*t(end);
    xo_con(i,:)=x(end,:);
    t2(i,:)=t(i,:);
       t_beg(i)=i*60-(60);
    
    % Computing oxygen consumption rates
    J = fluxes(x(end,:),param,ExpType,StateType);
    JO2_r(i) = (J(3)/2)*CF; %1.8957e5; % Unit conversion factor 
    end
 
elseif fig_flag==1
    
valy=ones(1,11) ;
valy=valy*30 ;
valx=[1;3;6;9;12;15;18;21;24;27;30;];
valx=[valx',valy]./1000; %conversion from mM to M
 
% Initialize 
  % Adding glutamate and malate
    xo(iGLU_c) = 0.01; %5e-3; 10 mM glutamate
    xo(iMAL_c) = 0.0025; %5e-3; 25 mM
    xo(iPi_c)  = 20e-3;
    xo(iADP_c) = 0.2e-3;  
    J_AtC = .015e-3; %M/sec/mg
    [t,x] = ode15s(@IsoMito_dXdT,[0 60],xo,options,param,ExpType,StateType,J_AtC);
    xo1 = x(end,:);
    
         
valy=ones(1,11) ;
valy=valy*30 ;
valx=[1;3;6;9;12;15;18;21;24;27;30;];
valx=[valx',valy]./1000;
    
    for i=1:1:11
     xo1(iPi_c)  = 20e-3;
     xo1(iADP_c) = 0.2e-3;   
    PCr_a(i) = valx(i);
    xo1(iPCr_c) = PCr_a(i); 
    J_AtC = .015e-3; %M/sec/mg
    [t,x] = ode15s(@IsoMito_dXdT,[0 60],xo1,options,param,ExpType,StateType,J_AtC);
    xo1 = x(end,:);

    H_x_PI(i) = x(end,iH_x);
    NADH_x_PI(i) = x(end,iNADH_x);
    Cred_i_PI(i) = x(end,iCred_i);
    dPsi_PI(i)  =x(end,idPsi);       
    ATP_x_PI(i) =x(end,iATP_x);
    ADP_x_PI(i) =x(end,iADP_x);
    ATP_c_PI(i) =x(end,iATP_c);
    ADP_c_PI(i) =x(end,iADP_c);
    PCr_c_PI(i) =x(end,iPCr_c);
    Pi_x_PI(i) =x(end,iPi_x);
    Pi_c_PI(i) =x(end,iPi_c);
  
    t_new(i) =i*t(end);
    t_beg(i)=i*t(1);
    xo_con(i,:)=x(end,:);
    t2(i,:)=t(i,:);
    
    % Computing oxygen consumption rates
    J = fluxes(x(end,:),param,ExpType,StateType);
    JO2_r(i) = (J(3)/2)*CF;%1.8957e5; % Unit conversion factor 
    end

%% Figure 4
elseif fig_flag ==3 
    
	valy=ones(1,11) ;
	valy=valy*30 ;
	valx=[1;3;6;9;12;15;18;21;24;27;30;];
	valx=[valx',valy]./1000; 
    
    % Adding glutamate and malate
    J_AtC = .015e-3; %M/sec/mg
    [t,x] = ode15s(@IsoMito_dXdT,[0 20],xo,options,param,ExpType,StateType,J_AtC);
    xo1=x(end,:);
    
    for i=1:1:5
    xo1(iPYR_c) = 0.005; %5e-3; 5 mM Pyruvate
    xo1(iMAL_c) = 0.0025; %5e-3; 2.5 mM
    PCr_a(i) = valx(i);
    xo1(iPCr_c) = PCr_a(i); 
    J_AtC = .015e-3; %M/sec/mg
    [t,x] = ode15s(@IsoMito_dXdT,[0 20],xo1,options,param,ExpType,StateType,J_AtC);
    
    t_new(i)=i*t(end);
    t_beg(i)=i*t(1);
    dPsi_PI(i) =x(end,idPsi);  
    Pi_c_PI(i) =x(end,iPi_c);
    ADP_c_PI(i)=x(end,iADP_c);
    ATP_c_PI(i) =x(end,iATP_c);
    PCr_c_PI(i) =x(end,iPCr_c);
    
    J = fluxes(x(end,:),param,ExpType,StateType);
    JO2_r(i) = (J(3)/2)*CF_mg;%1.8957e5; % Unit conversion factor 
    end
    
	%Observed data heart muscle fig. 4E    
	Obs_Jo2_PM_HM=[5729.70; 4549.05; 3682.12;3040.8;];
	Obs_dPsi_PM_HM=[152.52; 154.37; 156.08; 156.71;];
 

end    
    
%% Objective function to determine dG 
if (fig_flag == 1) || (fig_flag == 2)
load Table1
Y_obs= Table1.dG';
%[Oss,dG]=objfunctiondG(Y_obs,Pi_x_PI,ADP_x_PI,ATP_x_PI)

dG=[-33.8; -33.8762; -33.9929; -34.0881; -34.1951; -34.2850;...
    -34.3845; -34.4829; -34.5661; -34.6340; -34.7163;];
R= 8.3144598 ;
T= 273.13+37 ;
RT=R*T/1000  ;
Pi=Pi_c_PI;
    for i=1:1:11
dGATP_s(i) = (dG(i)) + RT * log((ADP_c_PI(i)*Pi(i))/ATP_c_PI(i));
    end
end

%% Loading state-2 experimental data curves (from Bose et al., 2003)
load Table1
%% Unit Conversion for plotting
ADP=(Table1.ADP);
ADP_c_PImm=ADP_c_PI*1000 ; %Conversion from M to mM
ATP_c_PImm=ATP_c_PI*1000 ;
PCr_c_PImm=PCr_c_PI*1000 ;


%% Plotting the results
if fig_flag == 1

figure(401); clf;
axes('position',[0.15 0.15 0.75 0.75]);
box on; hold on;
plot(Table1.dG,ADP,'ko','MarkerSize',12,'MarkerFaceColor',[1 1 1]);
plot(dGATP_s,ADP_c_PImm,'-k','LineWidth',2);
set(gca,'Fontsize',18);
xlabel('\Delta G_A_T_P (kJ/mol)');
ylabel('ADP (mM)');
print -f401 -r600 -depsc Muoio_ADP;

figure(402); clf;
axes('position',[0.15 0.15 0.75 0.75]);
box on; hold on;
plot(Table1.dG,Table1.ATP,'ko','MarkerSize',12,'MarkerFaceColor',[1 1 1]);
plot(dGATP_s,ATP_c_PImm,'-k','LineWidth',2);
set(gca,'Fontsize',18);
ylim([0 6]);
xlabel('\Delta G_A_T_P (kJ/mol)');
ylabel('ATP (mM)');
print -f402 -r600 -depsc Bose_MVO2;

figure(403); clf;
axes('position',[0.15 0.15 0.75 0.75]);
box on;hold on;
plot(Table1.dG,Table1.PCr,'ko','MarkerSize',12,'MarkerFaceColor',[1 1 1]);
plot(dGATP_s,PCr_c_PImm,'-k','LineWidth',2);set(gca,'Fontsize',18);
xlabel('\Delta G_A_T_P (kJ/mol)');
ylabel('PCr (mM)');


%Figure 1E 
Obs_Jo2_PM_HM=[1727; 2020;2313; 2636; 3165; 3899; 4898;6308;8365;10624;];
Obs_dG_PM_HM=[-14.9; -14.78; -14.68; -14.60; -14.497; -14.353; -14.19; -13.95; -13.65; -12.95;];
Obs_dG_PM_HM=Obs_dG_PM_HM*4.184;



elseif fig_flag ==2 
t_o2=[0.09;0.83; 2.48; 3.08; 3.94;5.35; 6.01; 8.34; 10.55; 12.65; 15.62;]; %minutes
t_o2=t_o2*60; % seconds
O2=[17.33; 19.07; 208.66; 218.63; 224.00; 225.67; 224.68; 181.53; 141.13; 107.13; 96.76;];

figure(405); clf;
axes('position',[0.15 0.15 0.75 0.75]);
box on;hold on;
plot(t_beg,JO2_r,'-k','LineWidth',2);
set(gca,'Fontsize',18);
plot(t_o2,O2,'ko','MarkerSize',12,'MarkerFaceColor', [1 1 1]);
xlabel('Time (sec)');
ylabel('JO_2 (pmol/sec/mL)');


end


%% End of the function



