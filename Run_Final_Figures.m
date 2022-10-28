%% Run Combined Model of Oxphos,TCA cycle and Beta-Oxdiation
%  This file runs simulations and parameter fitting for both the in vitro
%  and in vivo versions of the in silco mitchondrial metabolism model
%  representing cardiac muscle energetics

% Set default plot settings
set(groot,'defaultAxesFontSize',16); set(groot,'defaultLineLineWidth',2);

%% In vitro verification and parameter fitting

% Plot Figure 3: 
Reproduce_InVitro % Fig C-F
Plot_variability  % Fig 3A and B


saveas(1,'Fig3A.pdf');
saveas(2,'Fig3B.pdf');
saveas(402,'Fig3C.pdf'); 
saveas(403,'Fig3D.pdf'); 
saveas(401,'Fig3E.pdf');
saveas(405,'Fig3F.pdf');

clear all
close all

%% In Vivo Simulations
% Run setup Define Globals
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
	   i_J_AKGt i_J_SUCt i_J_FUMt i_J_GLUt i_J_ASPt i_J_CKe i_J_AKe i_J_AtC...
       iC16Carn_m iC16CoA_m iC16EnoylCoA_m iC16OHCoA_m iC16KetoCoA_m...   
       iC14Carn_cy iC14Carn_m iC14CoA_m iC14EnoylCoA_m iC14OHCoA_34e4r5m...    
       iC14KetoCoA_m iC12Carn_cy iC12Carn_m iC12CoA_m iC12EnoylCoA_m iC12OHCoA_m...   
       iC12KetoCoA_m iC10Carn_cy iC10Carn_m iC10CoA_m iC10EnoylCoA_m iC10OHCoA_m...  
       iC10KetoCoA_m iC8Carn_cy iC8Carn_m iC8CoA_m iC8EnoylCoA_m iC8OHCoA_m...  
       iC8KetoCoA_m iC6Carn_cy iC6Carn_m iC6CoA_m iC6EnoylCoA_m iC6OHCoA_m...
       iC6KetoCoA_m iC4Carn_cy iC4Carn_m iC4CoA_m...
       iC4EnoylCoA_m iC4OHCoA_m iC4KetoCoA_m iAcetylCoAMAT iFADH_m iNADHm iCoAMAT...
       iC16AcylCoACYT iC16Carn_cy;
% B-ox Vmax variables [vfcact; vrcact; vcpt2; vvlcad; vlcad; vmcad; vscad; vcrot; vmschad; vmckat; vmpt; vcpt1]
 
%% Setup Run
% Run setup file to set constants and initial conditions
setup; 
options=odeset('NonNegative', 1:62);
xo = xo_con;         % Initial Conditions
SF=10^6;             % Scaling Factor to Convert uM to M
iflag=100;           % iflag sets the ATP hydrolysis rate (work rate) to a fixed value
xo(109) = 2.68E-4;   % Palmitoyl-CoA initial cocndition: 2.68E-4*10 Fasted; % 26.8 uM Fed
xo(iPYR_c)= 120E-6;  % Pyruvate initial condition: 120E-6./10 Fasted; % At 120 uM 50% of ATP is supplied by Pyruvate

%% Figure 4A
% Plot Substrate Selection at rest; option to run after moderate exercise
% conditions and under fasted or fed conditions
     load('Varstruc_default');
     Varstruc.NADtot_p44=param(44);
    
     % Healthy: Set IC representative of 20 yr old (Wu et al., 2008)
     xo1 =xo; %
     xo1(iPI_c) = 5e-3;
     xo1(iATP_c)=0.011279;
     xo1(iPCr_c)=0.02198;
     xo1(iCr_c)=0.03917; 
     xo1(109) = 2.68E-5;  % Palmitoyl-CoA Concentration: 26.8 uM Fed  % 2.68E-4*10 Fasted;  
     xo1(iPYR_c)= 120E-5; % Pyruvate Concentration:  12  uM Fed       % 120E-6./10 Fasted;          

     Run_substrate_select(xo1,param,Varstruc,0.5e-3) % Rest
     %Run_substrate_select(xo1,param,Varstruc,1.0e-3)% Moderate Exercise
     colormap(bone);
     
     
     % Save figure to figure folder
     saveas(gcf,'Fig4A.pdf'); 
  
%% Figures 4B,4C,4D
     
     % Healthy: Set IC representative of 20 yr old (Wu et al., 2008)
     xo1 =xo; 
     xo1(iPI_c) = 5e-3;
     xo1(iATP_c)=0.011279;
     xo1(iPCr_c)=0.02198;
     xo1(iCr_c)=0.03917; 
     xo1(109) = 2.68E-4;  % Fasted State - Set Palmitoyl-CoA Concentration.  
     xo1(iPYR_c)= 120E-6; % Fasted State - Set Pyruvate Concentration
     param1=param; 
     
     sig=0.2;                                         % Add variability to sensitive parameters
     nm = 25;                                         % # of parameter sets in LHS
     Ind=[5; 9; 23; 35; 39; 12; 1; ];                 % Indices for parameters (VmaxAKG, VmaxMDH, Ki_NADHH_akg, ANT activity,Proton Leak Activity, CPT1 transporter, KmCpt1) 
     nInd = numel(Ind);
     ub= [param(Ind(1:nInd)) + param(Ind(1:nInd)).*sig, Varstruc.boxVmax(12)+ Varstruc.boxVmax(12).*sig, Varstruc.Kmcpt1CarCYT + Varstruc.Kmcpt1CarCYT.*sig ]; 
     lb= [param(Ind(1:nInd)) - param(Ind(1:nInd)).*sig, Varstruc.boxVmax(12)- Varstruc.boxVmax(12).*sig, Varstruc.Kmcpt1CarCYT - Varstruc.Kmcpt1CarCYT.*sig ];  
     % Generate LHS
     X_lhs=lhsdesign(nm,nInd+2);
     X_2= bsxfun(@plus,lb,bsxfun(@times,X_lhs,(ub-lb)));
          
     % Run simulation under fasted conditions
     [sim_pop]=Generate_Fig4b(Ind,X_2,xo1,param,Varstruc);
   
     % Run simulatin under fed conditions    
     % xo1(109) = 2.68E-6;  % Fed State - Set Palmitoyl-CoA Concentration.   
     % xo1(iPYR_c)= 1.2E-3; % Fed State - Set Pyruvate Concentration
     
   
     % [sim_pop2]=Generate_Fig4b(Ind,X_2,xo1,param,Varstruc);
     
     fig_num= 1; 
     Plot_Pred_Intv(sim_pop.pcr_atph,1,fig_num)   %fasted
     %Plot_Pred_Intv(sim_pop2.pcr_atph,2,fig_num) %fed
     ylim([0 2.5]);
     
     fig_num= 2; 
     Plot_Pred_Intv(sim_pop.pdhcsh,1,fig_num)   %fasted
     %Plot_Pred_Intv(sim_pop2.pdhcsh,2,fig_num) %fed
     ylim([0 1]);
    
     fig_num= 3; 
     Plot_Pred_Intv(sim_pop.pih,1,fig_num)   %fasted
     %Plot_Pred_Intv(sim_pop2.pih,2,fig_num) %fed
     ylim([0 6]);
     
     % Save figures to figure folder
     saveas(1,'Fig4B.pdf');
     saveas(2,'Fig4C.pdf');
     saveas(3,'Fig4D.pdf');
     
%% Simulating Mitochondrial Dysfunction Associated with Heart Failure
% Intervention 1: None - Healthy or Disease state 
% Intervention 2: DCA - increase PDH activity (fix aa to 1) 
% Intervention 3: Nicatinamide Riboside (NR) - Increase NAD by 150% (maximal effect)
% Intervention 4: Trimetazidine - inhibits vmckat  
%  Variable structure to pass through interventions or disease states
%  Varstruc.NADtot_p44=param(44);    % Total NAD pools
%  Varstruc.PDH_DR=0.20;             % Down-regulation of PDH activity
%  Varstruc.PDH_UR=1;                % Up-regulation of PDH activity
%  Varstruc.boxVmax=boxVmax1./10;    % Beta-oxidation Vmax parameters
%  Varstruc.PDH_aa=1;                % PDH activity
%  Varstruc.mito_rho=0.9;            % Density of mitochondria in the cell 
%  Varstruc.Trimetazidine=0;         % Concentration of Vmckat Inhibitor Trimetazidine

%%   Reduced Nucleotide Pools
     % Set IC representative of 80 yr old
     load('Varstruc_default');
     Varstruc.NADtot_p44=param(44);
     xo3 =xo;
     xo3(iPI_c) = 5e-3;
     xo3(iATP_c)= 0.002692368; 
     xo3(iPCr_c)=0.00948896;
     xo3(iCr_c)=0.009949086; 
     
     [sim_pop3]=Generate_Fig4b(Ind,X_2,xo3,param1,Varstruc);
      
     % Simulate Reduced Nucleotide Pools treated with dicholoroacetate
     % (DCA) Mechanism of action (MoA) is to increase flux through PDH.
     load('Varstruc_default');
     Varstruc.NADtot_p44=param(44);
     Varstruc.PDH_UR=1; % Increases flux through PDH (DCA MoA)
     xo3 =xo;
     xo3(iPI_c) = 5e-3;
     xo3(iATP_c)= 0.002692368; %0.0041235;
     xo3(iPCr_c)=0.00948896; %0.011571;
     xo3(iCr_c)=0.009949086; %0.01482;
     
     [sim_pop4]=Generate_Fig4b(Ind,X_2,xo3,param1,Varstruc);
      
      
%% Downregulate Beta oxidation and PDH activity

     % Set IC representative of 20 yr old and decrease PDH activity and
     % reduce Vmax of beta-oxidation enzymes to 1/10 
     load('Varstruc_default');
     Varstruc.NADtot_p44=param(44);
     Varstruc.PDH_DR=0.15;
     Varstruc.boxVmax=boxVmax./10;

     xo1 =xo; 
     xo1(iPI_c)  = 5e-3;
     xo1(iATP_c) = 0.011279;
     xo1(iPCr_c) = 0.02198;
     xo1(iCr_c)  = 0.03917; 
     
     [sim_pop5]=Generate_Fig4b(Ind,X_2,xo1,param1,Varstruc);
      
     % Simulate Decreased beta-oxidation and PDH activity treated with DCA
     load('Varstruc_default');
     Varstruc.NADtot_p44=param(44);
     Varstruc.PDH_UR = 1; 
     Varstruc.PDH_DR=0.15;
     Varstruc.boxVmax=boxVmax./10;
     
     param1 = param; 
     xo1 =xo; 
     xo1(iPI_c)  = 5e-3;
     xo1(iATP_c) = 0.011279;
     xo1(iPCr_c) = 0.02198;
     xo1(iCr_c)  = 0.03917; 
     
     [sim_pop6]=Generate_Fig4b(Ind,X_2,xo1,param1,Varstruc);
     
     
%% Decrease Mitochondrial Volume
    
     % Set IC representative of 20 yr old and decrease volume of mitochondria by 20%
     load('Varstruc_default');
     Varstruc.NADtot_p44	= param(44);
     Varstruc.mito_rho		= 0.80; % Set mito_rho to 80% of original value
  
     
     param1=param;
     xo1 =xo;
     xo1(iPI_c)   = 5e-3;
     xo1(iATP_c)  = 0.011279;
     xo1(iPCr_c)  = 0.02198;
     xo1(iCr_c)   = 0.03917;
     
    [sim_pop7]=Generate_Fig4b(Ind,X_2,xo1,param1,Varstruc);
 
    
     % Simulate decrease volume of mitochondria by 20% treated with DCA
     load('Varstruc_default');
     Varstruc.NADtot_p44	= param(44);
      Varstruc.mito_rho		= 0.80;	       % Set mito_rho to 80
     Varstruc.PDH_UR 		= 1;           % MoA of DCA
 
     param1 = param; 
     xo1 =xo; 
     xo1(iPI_c)   = 5e-3;
     xo1(iATP_c)  = 0.011279;
     xo1(iPCr_c)  = 0.02198;
     xo1(iCr_c)   = 0.03917; 
     
     [sim_pop8]=Generate_Fig4b(Ind,X_2,xo1,param1,Varstruc);
     

%% Plot Figure 5
    
% 0.5 mM/sec ATP hydrolysis
y(1)=mean(sim_pop.pcr_atph([2:end],2)); 
y(2)=mean(sim_pop5.pcr_atph([2:end],2));
y(3)=mean(sim_pop7.pcr_atph([2:end],2));
y(4)=mean(sim_pop3.pcr_atph([2:end],2));

% 1.0 mM/sec ATP hydrolysis
m(1)=mean(sim_pop.pcr_atph([2:end],7));
m(2)=mean(sim_pop5.pcr_atph([2:end],7));
m(3)=mean(sim_pop7.pcr_atph([2:end],7));
m(4)=mean(sim_pop3.pcr_atph([2:end],7));

vv=[sim_pop.pcr_atph(2:end,2),sim_pop5.pcr_atph(2:end,2),sim_pop7.pcr_atph(2:end,2),sim_pop3.pcr_atph(2:end,2)];
vx=[sim_pop.pcr_atph(2:end,7),sim_pop5.pcr_atph(2:end,7),sim_pop7.pcr_atph(2:end,7),sim_pop3.pcr_atph(2:end,7)];




z=[y;m];
vv1=vv;%./vv(1,:);
vx1=vx;%./vx(1,:);

% Calculate 95% prediction interval
for k=1:4
    ub(1,k)=prctile(vv1(:,k),95)-mean(vv1(:,k));
    lb(1,k)=-prctile(vv1(:,k),5)+mean(vv1(:,k));
    ub(2,k)=prctile(vx1(:,k),95)-mean(vx1(:,k));
    lb(2,k)=-prctile(vx1(:,k),5)+mean(vx1(:,k));
end

name ={'0.5 mM/s', '1.0 mM/s'};

figure;
ym=bar(z,'grouped','FaceColor',"flat");
for k = 1:size(ym,2)
    ym(k).CData = k;
end
colormap(bone); 

set(gca,'xticklabel',name); 
xlabel('ATP Hydrolysis Rate');ylabel({' PCr/ATP Ratio'; '(fraction of control)'});

ngroups = size(z, 1);
nbars = size(z, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for j = 1:1:nbars
    % Calculate center of each bar
xpos = (1:ngroups) - (groupwidth/2) + (2*j-1) * (groupwidth / (2*nbars));
hold on
    errorbar(xpos, z(:,j), lb(:,j),ub(:,j), 'k', 'linestyle', 'none');
hold off
end

lgd=legend('Control','Case 1: \downarrow PDH \beta-Ox','Case 2: \downarrow V_m_i_t_o',...
'Case 3: \downarrow Nucleotides','90% Prediction Interval');
lgd.Location = 'northeast'; lgd.FontSize=12;
ylim([0  2.5]); yticks([0  0.5  1 1.5 2 2.5]);

saveas(gcf,'Fig5.pdf');
%% Plot Figure 6
Np=linspace(0.36,1.3,10);
figure;
plot(Np,mean(sim_pop.pih([2:end],:)),'k-')
hold on
mean_pop5=mean(sim_pop5.pih([2:7],:));
plot(Np(2:7),mean_pop5(2:7),'k-.')
plot(Np,mean(sim_pop6.pih([2:end],:)),'k:')

xx1=[0.3 1.3]; yy1=[2 2]; plot(xx1,yy1,'--','color',[0.5 0.5 0.5]);  % Line denotes fatigue
xlabel('ATP Hydrolysis (mM/s)'); ylabel('Inorganic Phosphate (mM)');
title('Case 1'); xlim([0.2 1.4]); ylim([0  10]);
lg=legend('Healthy','\downarrow PDH \beta-Ox','\downarrow PDH \beta-Ox + DCA','Fatigue (2mM Pi)');
lg.Location = 'northwest';

saveas(gcf,'Fig6A.pdf');

%% Save Data

simpop.sim1=sim_pop;    % Healthy fasted state
%simpop.sim2=sim_pop2;   % Healthy fed state
simpop.sim3=sim_pop3;   % Reduced nucleotides
simpop.sim4=sim_pop4;   % Reduced nucleotides + DCA
simpop.sim5=sim_pop5;   % Down-regulate PDH and beta-oxidation
simpop.sim6=sim_pop6;   % Down-regulate PDH and beta-oxidation + DCA
simpop.sim7=sim_pop7;   % Reduced mitochondria volume
simpop.sim8=sim_pop8;   % Reduced mitochondria volume + DCA

save('simpop');

% Figure 6B
Fig6B();
    
